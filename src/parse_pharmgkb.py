#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Process PAnno data.
"""

# Import modules
import os, json, html, toml, re, itertools
import pandas as pd
import numpy as np
from utils.pgi_functions import connect_database, pymysql_cursor, load_logging_cfg
from utils.parse_clinann import infer_response

cfile = './conf/panno.toml'
cf = toml.load(cfile, _dict=dict)
DEFAULT = cf['DEFAULT']
RAW = cf['RAW']
TSV = cf['TSV']

db, cursor = connect_database(DEFAULT)
logger = load_logging_cfg(DEFAULT)


### Some gadgets
def columns_length(df):
  for i in range(0, df.shape[1]):
    print(df.columns.to_list()[i])
    print(df.iloc[:,i].astype(str).str.len().max())



#########################################
## Primary Data
#########################################
##### PharmGKB.Gene
genes = pd.read_csv(RAW["genes"], sep="\t", dtype="object")
#ann_genes = genes[(genes["Is VIP"] == "Yes") | (genes["Has Variant Annotation"] == "Yes") | (genes["Has CPIC Dosing Guideline"] == "Yes")]
Gene = genes[['Symbol', 'Name', 'PharmGKB Accession Id', 'NCBI Gene ID', 'HGNC ID', 'Ensembl Id', 'Chromosome', 'Chromosomal Start - GRCh37', 'Chromosomal Stop - GRCh37', 'Chromosomal Start - GRCh38', 'Chromosomal Stop - GRCh38']].reset_index(drop = True)

table_p = TSV["Gene"]
Gene.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Gene FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Symbol, Name, PAID, EntrezID, HGNCID, EnsemblID, Chr, GRCh37Start, GRCh37End, GRCh38Start, GRCh38End);" % table_p)
row_count

# Insert the '-'
pymysql_cursor('INSERT INTO Gene(Symbol) VALUES("-");')


##### PharmGKB.Variant
variants = pd.read_csv(RAW["variants"], sep="\t", dtype="object")
Variant = variants[['Variant Name', 'Location', 'Variant ID']]
Variant["Type"] = "Variant"

# Add Haplotype into PharmGKB.Variant table as a specific variant
var_drug_ann = pd.read_csv(RAW["var_drug_ann"], sep="\t", usecols=["Variant/Haplotypes", "Gene"])
var_fa_ann = pd.read_csv(RAW["var_fa_ann"], sep="\t", usecols=["Variant/Haplotypes", "Gene"])
var_pheno_ann = pd.read_csv(RAW["var_pheno_ann"], sep="\t", usecols=["Variant/Haplotypes", "Gene"])
clinical_annotations = pd.read_csv(RAW["clinical_annotations"], sep="\t", usecols=["Variant/Haplotypes", "Gene"])

merge = pd.concat([var_drug_ann, var_fa_ann, var_pheno_ann, clinical_annotations], axis=0).drop_duplicates().reset_index(drop = True)
haplotype_merge = merge[merge["Variant/Haplotypes"].str.extract(r'(rs\d+)', expand = False).isnull()]
Haplotype = haplotype_merge.drop(['Variant/Haplotypes'], axis=1).join(haplotype_merge['Variant/Haplotypes'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('Variant Name')).drop_duplicates().reset_index(drop = True)
Haplotype = pd.concat([Haplotype, pd.DataFrame(columns = ['Location', 'Variant ID'])])
Haplotype["Type"] = "Haplotype"

# Merge Variant and Haplotype
VariantOrHaplotype = pd.concat([Variant, Haplotype[['Variant Name', 'Location', 'Variant ID', 'Type']]], axis=0)

table_p = TSV["Variant"]
VariantOrHaplotype.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Variant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Name, NCID, PAID, Type);" % table_p)
row_count


##### PharmGKB.VariantSynonyms
variants = pd.read_csv(RAW["variants"], sep="\t")
VariantSynonyms = variants[['Variant ID', 'Synonyms']].reset_index(drop = True)
for index, row in VariantSynonyms.iterrows():
  VariantSynonyms.loc[index, "Variant ID"] = pymysql_cursor("SELECT ID FROM Variant Where PAID = '%s';" % row["Variant ID"])

VariantSynonyms = VariantSynonyms.drop(['Synonyms'], axis=1).join(VariantSynonyms['Synonyms'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('Synonym'))
columns_length(VariantSynonyms)

table_p = TSV["VariantSynonyms"]
VariantSynonyms.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE VariantSynonyms FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (VariantID, Synonym);" % table_p)
row_count


##### PharmGKB.GeneHasVariant
variants = pd.read_csv(RAW["variants"], sep="\t")
GeneHasVariant = variants[variants["Gene IDs"].isnull() == False].copy()
GeneHasVariant = GeneHasVariant[["Gene IDs", "Variant ID"]].reset_index(drop = True)
GeneHasVariant = GeneHasVariant.drop(['Gene IDs'], axis=1).join(GeneHasVariant['Gene IDs'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('GeneID'))
GeneHasVariant = GeneHasVariant.rename(columns = {"Variant ID": "VariantID"})
for index, row in GeneHasVariant.iterrows():
  GeneHasVariant.loc[index, "GeneID"] = pymysql_cursor("SELECT ID FROM Gene Where PAID = '%s';" % row["GeneID"])
  GeneHasVariant.loc[index, "VariantID"] = pymysql_cursor("SELECT ID FROM Variant Where PAID = '%s';" % row["VariantID"])

# Add GeneHasHaplotype
GeneHasHaplotype = Haplotype[["Variant Name", "Gene"]].copy()
for index, row in GeneHasHaplotype.iterrows():
  GeneHasHaplotype.loc[index, "Variant Name"] = pymysql_cursor("SELECT ID FROM Variant Where Name = '%s';" % row["Variant Name"])
  GeneHasHaplotype.loc[index, "Gene"] = pymysql_cursor("SELECT ID FROM Gene Where Symbol = '%s';" % row["Gene"])

GeneHasHaplotype = GeneHasHaplotype.rename(columns = {"Gene": "GeneID", "Variant Name": "VariantID"})
GeneHasVariantOrHaplotype = pd.concat([GeneHasVariant, GeneHasHaplotype], axis=0)

table_p = TSV["GeneHasVariant"]
GeneHasVariantOrHaplotype.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE GeneHasVariant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (VariantID, GeneID);" % table_p)
row_count


##### PharmGKB.Drug
drugs = pd.read_csv(RAW["drugs"], sep="\t")
chemicals = pd.read_csv(RAW["chemicals"], sep="\t")
Drug = pd.concat([drugs, chemicals], axis=0).drop_duplicates().reset_index(drop = True)
Drug.insert(0, "ID", (pd.DataFrame(Drug.index)+1)[0].to_list())

# Prepare for the next table
Drug[["Name", "Generic Names", "Trade Names"]] = Drug[["Name", "Generic Names", "Trade Names"]].replace(np.nan, "")
Drug["Synonyms"] = Drug["Name"].str.cat(Drug["Generic Names"], sep = ',"').str.cat(Drug["Trade Names"], sep = ',"')
DrugSynonyms = Drug[["ID", "Synonyms"]].rename(columns = {"ID": "DrugID"})

# Current table
Drug = Drug[["ID", "Name", "PharmGKB Accession Id", "Generic Names", "Trade Names", "Cross-references", "Type"]]

for index, row in Drug.iterrows():
  if row["Generic Names"] is np.nan:
    row["Generic Names"] = ""
  else:
    row["Generic Names"] = row["Generic Names"].split(',"')[0]
  if row["Trade Names"] is np.nan:
    row["Trade Names"] = ""
  else:
    row["Trade Names"] = row["Trade Names"].split(',"')[0]
  Drug.iloc[index] = row

#columns_length(Drug)
#Drug[Drug["PharmGKB Accession Id"].duplicated()]
table_p = TSV["Drug"]
Drug.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Drug FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (ID, Name, PAID, GenericName, TradeName, CrossRef, Type);" % table_p)
row_count


##### PharmGKB.DrugSynonyms
DrugSynonyms = DrugSynonyms.dropna().reset_index(drop = True)
DrugSynonyms = DrugSynonyms.drop(['Synonyms'], axis=1).join(DrugSynonyms['Synonyms'].str.split(',"', expand=True).stack().reset_index(level=1, drop=True).rename('Synonym').str.strip('"'))
DrugSynonyms = DrugSynonyms.drop_duplicates().reset_index(drop = True)

DrugSynonyms = DrugSynonyms[DrugSynonyms["Synonym"] != ""]
columns_length(DrugSynonyms)

table_p = TSV["DrugSynonyms"]
DrugSynonyms.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE DrugSynonyms FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (DrugID, Synonym);" % table_p)
row_count


##### PharmGKB.Phenotype
phenotypes = pd.read_csv(RAW["phenotypes"], sep="\t")
Phenotype = phenotypes[["PharmGKB Accession Id", "Name", "Alternate Names"]].drop_duplicates()
Phenotype.insert(0, 'ID', (pd.DataFrame(Phenotype.index)+1)[0].to_list())

table_p = TSV["Phenotype"]
Phenotype[["ID", "PharmGKB Accession Id", "Name"]].to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Phenotype FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (ID, PAID, Name);" % table_p)
row_count

pymysql_cursor('INSERT INTO Phenotype(Name, PAID, Type) VALUES("-", "-", "-");')


##### PharmGKB.PhenotypeSynonyms
Phenotype[["Name", "Alternate Names"]] = Phenotype[["Name", "Alternate Names"]].replace(np.nan, "")
Phenotype["Alternate Names"] = Phenotype["Name"].str.cat(Phenotype["Alternate Names"], sep = ',"')

PhenotypeSynonyms = Phenotype[["ID", "Alternate Names"]].drop_duplicates().rename(columns = {'ID': 'PhenotypeID'}).reset_index(drop = True)
PhenotypeSynonyms = PhenotypeSynonyms.drop(['Alternate Names'], axis=1).join(PhenotypeSynonyms['Alternate Names'].str.split(',"', expand=True).stack().reset_index(level=1, drop=True).rename('Synonym').str.strip('"'))
PhenotypeSynonyms = PhenotypeSynonyms.drop_duplicates()
PhenotypeSynonyms = PhenotypeSynonyms[PhenotypeSynonyms["Synonym"] != ""]
columns_length(PhenotypeSynonyms)

table_p = TSV["PhenotypeSynonyms"]
PhenotypeSynonyms.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PhenotypeSynonyms FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (PhenotypeID, Synonym);" % table_p)
row_count

# Determine if this phenotype is cancer
cancer_reg = "|".join(["cancer", "cancers", "tumors", "tumor", "tumours", "tumour", "carcinomas", "carcinoma", "neoplasms", "neoplasm", "epitheliomas", "epithelioma", "malignancies", "malignancy", "adenocarcinoma", "cystadenoma", "cystadenocarcinoma", "sarcoma"])
for index, row in PhenotypeSynonyms.iterrows():
  if re.findall(cancer_reg, row.Synonym.lower()):
    pymysql_cursor("UPDATE PharmGKB.Phenotype SET Type = 'Cancer' WHERE ID = '%s';" % row.PhenotypeID)
  else:
    pymysql_cursor("UPDATE PharmGKB.Phenotype SET Type = 'Not cancer' WHERE ID = '%s';" % row.PhenotypeID)



#########################################
## Annotation Related Dictionary
#########################################
##### PharmGKB.EvidenceTypeDic
clinical_ann_evidence = pd.read_csv(RAW["clinical_ann_evidence"], sep="\t")
EvidenceTypeDic = pd.DataFrame({"Name": list(set(clinical_ann_evidence["Evidence Type"].to_list()))})

table_p = TSV["EvidenceTypeDic"]
EvidenceTypeDic.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE EvidenceTypeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Name);" % table_p)
row_count


##### PharmGKB.PhenotypeCategoryDic
clinical_annotations = pd.read_csv(RAW["clinical_annotations"], sep="\t")
pheno_category = list(set(clinical_annotations["Phenotype Category"].str.split(";", expand=True).stack().to_list())) + ["PD"]
PhenotypeCategoryDic = pd.DataFrame({"Name": pheno_category})

table_p = TSV["PhenotypeCategoryDic"]
PhenotypeCategoryDic.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PhenotypeCategoryDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Name);" % table_p)
row_count



#########################################
## Variant Annotation
#########################################
##### PharmGKB.VarAnn
# var_drug_ann.tsv: CYP2C19􏰀17 -> CYP2C19*17
# 􏰃1.569 to <2.075 mg/kg/ day -> ≥1.569 to ≤2.075 mg/kg/day
# 􏰃1.54 to <2.05 mg/kg/day -> ≥1.54 to ≤2.05 mg/kg/day
# "The 15-year -> The 15-year
var_drug_ann = pd.read_csv(RAW["var_drug_ann"], sep="\t")
var_fa_ann = pd.read_csv(RAW["var_fa_ann"], sep="\t")
# var_pheno_ann.tsv: "Although roughly -> Although roughly
var_pheno_ann = pd.read_csv(RAW["var_pheno_ann"], sep="\t")

var_drug_ann["Type"] = "Variant Drug Annotation"
var_fa_ann["Type"] = "Variant Functional Assay Annotation"
var_pheno_ann["Type"] = "Variant Phenotype Annotation"

var_ann = pd.concat([var_drug_ann, var_fa_ann, var_pheno_ann], axis=0).drop_duplicates().reset_index(drop = True)
VarAnn = var_ann.drop(['Phenotype Category'], axis=1).join(var_ann['Phenotype Category'].str.split(',"', expand=True).stack().reset_index(level=1, drop=True).rename('Phenotype Category').str.strip('"')).reset_index(drop = True)

# Replace the nan into Other
nan_index = VarAnn[VarAnn['Phenotype Category'].isna()].index.to_list()
VarAnn.loc[nan_index, 'Phenotype Category'] = "Other"

for index, row in VarAnn.iterrows():
  row["Type"] = pymysql_cursor("SELECT ID FROM EvidenceTypeDic WHERE Name ='%s';" % row["Type"])
  row["Phenotype Category"] = pymysql_cursor("SELECT ID FROM PhenotypeCategoryDic WHERE Name ='%s';" % row["Phenotype Category"])
  VarAnn.iloc[index] = row

table_p = TSV["VarAnn"]
VarAnn.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor('LOAD DATA LOCAL INFILE "%s" INTO TABLE VarAnn FIELDS TERMINATED BY "\\t" LINES TERMINATED BY "\\n" IGNORE 1 LINES (AnnotationID, VariantOrHaplotype, Gene, Drug, PMID, Significance, Note, Sentence, Allele, Specialty, TypeID, CategoryID);' % table_p)
row_count




#########################################
## Clinical Annotation
#########################################
##### PharmGKB.ClinAnnEvidence
clinical_ann_evidence = pd.read_csv(RAW["clinical_ann_evidence"], sep="\t", dtype="object")
for index, row in clinical_ann_evidence.iterrows():
  row["Evidence Type"] = pymysql_cursor("SELECT ID FROM EvidenceTypeDic WHERE Name ='%s';" % row["Evidence Type"])
  clinical_ann_evidence.iloc[index] = row

clinical_ann_evidence["PMID"] = clinical_ann_evidence["PMID"].replace(np.nan, '')
table_p = TSV["ClinAnnEvidence"]
clinical_ann_evidence.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnnEvidence FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (CAID, EvidenceID, TypeID, URL, PMID, Summary, Score);" % table_p)
row_count


##### PharmGKB.ClinAnn
clinical_annotations = pd.read_csv(RAW["clinical_annotations"], sep="\t")
clinical_ann_alleles = pd.read_csv(RAW["clinical_ann_alleles"], sep="\t")
# Merge annotations and alleles
clinical_ann = pd.merge(clinical_annotations, clinical_ann_alleles, on="Clinical Annotation ID")
clinical_ann.columns
# Remove level 4
clinical_ann = clinical_ann[clinical_ann['Level of Evidence'] != '4'].reset_index(drop = True)

# Replace the nan into -
nan_index = clinical_ann[clinical_ann['Gene'].isna()].index.to_list()
clinical_ann.loc[nan_index, 'Gene'] = "-"
nan_index = clinical_ann[clinical_ann['Phenotype Category'].isna()].index.to_list()
clinical_ann.loc[nan_index,'Phenotype Category'] = "Other"
# Add columns
clinical_ann = pd.concat([clinical_ann, pd.DataFrame(columns = ['Allele1', 'Allele2', 'Annotation1', 'Annotation2', 'Function1', 'Function2', 'Score1', 'Score2', 'CPICPhenotype', 'PAnnoPhenotype'])], axis=1)


# Splite the variants and phenotype categories
# clinical_ann = clinical_ann.drop(['Variant/Haplotypes'], axis=1).join(clinical_ann['Variant/Haplotypes'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('Variant/Haplotypes')).reset_index(drop = True)
clinical_ann = clinical_ann.drop(['Gene'], axis=1).join(clinical_ann['Gene'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Gene')).reset_index(drop = True)
clinical_ann = clinical_ann.drop(['Phenotype Category'], axis=1).join(clinical_ann['Phenotype Category'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Phenotype Category')).reset_index(drop = True)

### ! Don't need to remove duplications after don't split Variant/Haplotypes
# Remove some duplications happened when mergeing.
# For example, clinical_ann[clinical_ann['Clinical Annotation ID'] == 1184746746]
# retained_index = []
for index, row in clinical_ann.iterrows():
  gene = row['Gene']
  alleles = row['Genotype/Allele']
  # Cleaning data
  if type(alleles) == str:
    alleles = re.sub('A-968C_376G', 'A- 968C_376G', re.sub('A-680T_376G', 'A- 680T_376G', re.sub('A-202A_376G', 'A- 202A_376G', alleles)))
    clinical_ann.loc[index, 'Genotype/Allele'] = alleles
  # haplotype = row['Variant/Haplotypes'].replace(gene, '').strip()
  # # Single locus
  # if haplotype.startswith('rs') or type(alleles) != str:
  #   retained_index.append(index)
  # # Others
  # elif haplotype in alleles:
  #   retained_index.append(index)

# clinical_ann.shape; len(retained_index) #(30459, 26), (17001, 26)
# clinical_ann = clinical_ann.iloc[retained_index,].reset_index(drop = True)
clinical_ann.shape # (16324, 26)
clinical_ann = clinical_ann.drop(['Latest History Date (YYYY-MM-DD)', 'URL'], axis=1)


################# Split clinical_ann into different classes #################
table_cols = ['CAID', 'Gene', 'Variant', 'Allele1', 'Allele2', 'Annotation1', 'Annotation2', 'Function1', 'Function2', 'Score1', 'Score2', 'CPICPhenotype', 'PAnnoPhenotype', 'Drug', 'Phenotypes', 'EvidenceLevel', 'LevelOverride', 'LevelModifier', 'Score', 'PMIDCount', 'EvidenceCount', 'Specialty', 'PhenotypeCategory']

## ! PAnno, Splite the data into THREE parts:
# ! 1. SNP/Indel: variants with dbSNP ID, genotypes only with clinical functions
# ! 2. diplotypes: drug-metabolizing enzymes with identified functions
# ! 3. HLA related genes
## Class 1: SNP/Indel
# Some star alleles in CPIC may be presented in the form of rsID in PharmGKB
# Therefore, searching the AlleleManual to find the CPIC phenotype
small_variants = clinical_ann[clinical_ann['Variant/Haplotypes'].str.startswith('rs')].reset_index(drop = True)
for index, row in small_variants.iloc[11102:13379,].iterrows():
  # Splite 'Genotype/Allele' into Allele1 and Allele2
  genotype = row['Genotype/Allele']
  if len(genotype) == 2:
    alleles = [genotype[0], genotype[1]]
  else:
    alleles = genotype.split("/")
  # Allele1
  row.Allele1 = alleles[0]
  tmp1 = cursor.execute('SELECT Allele, Function, FunctionManual FROM AlleleFunctionality WHERE Gene = "%s" AND Variant = "%s" AND AlleleManual = "%s";' % (row.Gene, row['Variant/Haplotypes'], row.Allele1))
  tmp1 = cursor.fetchall()
  if tmp1:
    row.Function1 = tmp1[0][2]
    A1 = tmp1[0][0]
  # Allele2
  if len(alleles) == 1:
    row.Allele2 = ''
    A2 = ''
  else:
    row.Allele2 = alleles[1]
    tmp2 = cursor.execute('SELECT Allele, Function, FunctionManual FROM AlleleFunctionality WHERE Gene = "%s" AND Variant = "%s" AND AlleleManual = "%s";' % (row.Gene, row['Variant/Haplotypes'], row.Allele2))
    tmp2 = cursor.fetchall()
    if tmp2:
      row.Function2 = tmp2[0][2]
      A2 = tmp2[0][0]
  # Annotation
  row.Annotation1 = row['Annotation Text']
  row.Annotation2 = ''
  row.PAnnoPhenotype = infer_response(row['Annotation Text'], row['Phenotype Category'])
  response_score = {'Unknown': np.nan, 'Uncertain': np.nan, 'No': 0, 'Decreased': 0.5, 'Normal': 1, 'Increased': 2}
  row.Score1 = response_score[row.PAnnoPhenotype]
  # CPICPhenotype
  row.CPICPhenotype = pymysql_cursor('SELECT Phenotype FROM DiplotypePhenotype WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s";' % (row.Gene, A1, A2))
  # Update
  small_variants.iloc[index] = row

small_variants = small_variants.drop(['Genotype/Allele', 'Annotation Text', 'Allele Function'], axis=1)
small_variants = small_variants.rename(columns={'Clinical Annotation ID':'CAID', 'Variant/Haplotypes':'Variant', 'Level of Evidence':'EvidenceLevel', 'Level Override':'LevelOverride', 'Level Modifiers':'LevelModifier', 'PMID Count':'PMIDCount', 'Evidence Count':'EvidenceCount', 'Drug(s)':'Drug', 'Phenotype(s)':'Phenotypes', 'Specialty Population':'Specialty', 'Phenotype Category':'PhenotypeCategory'})
small_variants = small_variants[table_cols]
small_variants.to_csv('./data/pharmgkb/intermediate/ClinAnn_small_variants.txt', sep='\t', index=0)


## Class 2 & 3: diplotypes and HLA related genes
# Need to combinate by ourselves
diplotypes = clinical_ann[clinical_ann['Variant/Haplotypes'].str.startswith('rs') == False]
diplotypes = diplotypes.drop(['Genotype/Allele'], axis=1).join(diplotypes['Genotype/Allele'].str.split(' \+ ', expand=True).stack().reset_index(level=1, drop=True).rename('Genotype/Allele')).reset_index(drop = True)
# Remove GSTT1 non-null, GSTT1 null
diplotypes = diplotypes[diplotypes['Variant/Haplotypes'] != 'GSTT1 non-null, GSTT1 null']
# diplotypes[diplotypes['Annotation Text'].str.contains('Patients with osteosarcoma and a GSTT1 null genotype')]['Genotype/Allele'] = 'null'

D1 = diplotypes[diplotypes['Genotype/Allele'].str.contains('/') == False]
caids = D1['Clinical Annotation ID'].drop_duplicates().to_list()
diplotypes_new = pd.DataFrame(columns = diplotypes.columns); index = 0
for id in caids:
  sub_df = D1[D1['Clinical Annotation ID'] == id]
  phes = sub_df['Phenotype Category'].drop_duplicates().to_list()
  all_hap_raw = sub_df['Genotype/Allele'].drop_duplicates().to_list()
  all_hap = []
  for hap in all_hap_raw:
    if '/' not in hap:
      all_hap.append(hap)
  all_hap = list(set(all_hap))
  # Gene info
  Gene = sub_df['Gene'].drop_duplicates().to_list()[0]
  chr = pymysql_cursor('SELECT Chr FROM Gene WHERE Symbol = "%s";')
  for phe in phes:
    phe_df = sub_df[sub_df['Phenotype Category'] == phe]
    combs = list(itertools.combinations_with_replacement(all_hap, 2))
    for dip in combs:
      row = phe_df.iloc[0,:].copy()
      # Allele1
      row.Allele1 = dip[0]
      row.Annotation1 = phe_df[phe_df['Genotype/Allele'] == row.Allele1]['Annotation Text'].to_list()[0]
      row.Function1 = pymysql_cursor('SELECT FunctionManual FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (Gene, row.Allele1.replace('xN', 'x2')))
      infer_func1 = infer_response(row.Annotation1, phe)
      if row.Function1 is None:
        # row.Function1 = phe_df[phe_df['Genotype/Allele'] == row.Allele1]['Allele Function'].to_list()[0]
        row.Function1 = infer_func1
      # Allele2
      row.Allele2 = dip[1]
      row.Annotation2 = phe_df[phe_df['Genotype/Allele'] == row.Allele2]['Annotation Text'].to_list()[0]
      row.Function2 = pymysql_cursor('SELECT FunctionManual FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (Gene, row.Allele2.replace('xN', 'x2')))
      infer_func2 = infer_response(row.Annotation2, phe)
      if row.Function2 is None:
        # row.Function2 = phe_df[phe_df['Genotype/Allele'] == row.Allele2]['Allele Function'].to_list()[0]
        row.Function2 = infer_func2
      # CPICPhenotype
      row.CPICPhenotype = pymysql_cursor('SELECT Phenotype FROM DiplotypePhenotype WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s";' % (row.Gene, row.Allele1.replace('xN', 'x2'), row.Allele2.replace('xN', 'x2')))
      if row.CPICPhenotype is None:
        row.CPICPhenotype = pymysql_cursor('SELECT Phenotype FROM DiplotypePhenotype WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s";' % (row.Gene, row.Allele2.replace('xN', 'x2'), row.Allele1.replace('xN', 'x2')))
      # PAnnoPhenotype: calculating based on inferred functions
      response_score = {'Unknown': np.nan, 'Uncertain': np.nan, 'No': 0, 'Decreased': 0.5, 'Normal': 1, 'Increased': 2}
      row.Score1 = response_score[infer_func1]; row.Score2 = response_score[infer_func2]
      total_score = row.Score1 + row.Score2
      if total_score in [0, 0.5, 1, 1.5]:
        row.PAnnoPhenotype = 'Decreased'
      elif total_score in [2, 2.5]:
        row.PAnnoPhenotype = 'Normal'
      elif total_score in [3, 4]:
        row.PAnnoPhenotype = 'Increased'
      else:
        row.PAnnoPhenotype = 'Indeterminate'
      # Update
      diplotypes_new.at[index, :] = row
      index = index+1


D2 = diplotypes[diplotypes['Genotype/Allele'].str.contains('/')].reset_index(drop = True)
for index, row in D2.iterrows():
  Gene = row.Gene
  dip = row['Genotype/Allele'].split("/")
  # Allele1
  row.Allele1 = dip[0]
  row.Annotation1 = row['Annotation Text']
  row.Function1 = pymysql_cursor('SELECT FunctionManual FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (Gene, row.Allele1.replace('xN', 'x2')))
  infer_func1 = infer_response(row.Annotation1, phe)
  if row.Function1 is None:
    row.Function1 = infer_func1
  # Allele2
  row.Allele2 = dip[1]
  row.Annotation2 = ''
  row.Function2 = pymysql_cursor('SELECT FunctionManual FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (Gene, row.Allele2.replace('xN', 'x2')))
  infer_func2 = infer_response(row.Annotation2, phe)
  if row.Function2 is None:
    row.Function2 = infer_func2
  # CPICPhenotype
  row.CPICPhenotype = pymysql_cursor('SELECT Phenotype FROM DiplotypePhenotype WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s";' % (row.Gene, row.Allele1.replace('xN', 'x2'), row.Allele2.replace('xN', 'x2')))
  if row.CPICPhenotype is None:
    row.CPICPhenotype = pymysql_cursor('SELECT Phenotype FROM DiplotypePhenotype WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s";' % (row.Gene, row.Allele2.replace('xN', 'x2'), row.Allele1.replace('xN', 'x2')))
  # PAnnoPhenotype
  row.PAnnoPhenotype = infer_response(row.Annotation1, phe)
  row.Score1 = response_score[infer_func1]; row.Score2 = response_score[infer_func2]
  D2.iloc[index] = row

diplotypes_merge = pd.concat([diplotypes_new, D2], axis = 0)
diplotypes_out = diplotypes_merge.drop(['Genotype/Allele', 'Annotation Text', 'Allele Function'], axis=1)
diplotypes_out = diplotypes_out.rename(columns={'Clinical Annotation ID':'CAID', 'Variant/Haplotypes':'Variant', 'Level of Evidence':'EvidenceLevel', 'Level Override':'LevelOverride', 'Level Modifiers':'LevelModifier', 'PMID Count':'PMIDCount', 'Evidence Count':'EvidenceCount', 'Drug(s)':'Drug', 'Phenotype(s)':'Phenotypes', 'Specialty Population':'Specialty', 'Phenotype Category':'PhenotypeCategory'})
diplotypes_out = diplotypes_out[table_cols]
diplotypes_out.to_csv('./data/pharmgkb/intermediate/ClinAnn_diplotypes.txt', sep='\t', index=0)


# Merge the above dataframe
clin_ann_snpindel = pd.read_csv('./data/pharmgkb/intermediate/ClinAnn_small_variants.txt', sep='\t')
clin_ann_diplotype = pd.read_csv('./data/pharmgkb/intermediate/ClinAnn_diplotypes.txt', sep='\t')
clin_ann_merge = pd.concat([clin_ann_snpindel, clin_ann_diplotype], axis=0).reset_index(drop = True)

# Splite the drugs
clin_ann_merge = clin_ann_merge.drop(['Drug'], axis=1).join(clin_ann_merge['Drug'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Drug')).reset_index(drop = True)

table_cols = ['CAID', 'Gene', 'Variant', 'Allele1', 'Allele2', 'Annotation1', 'Annotation2', 'Function1', 'Function2', 'Score1', 'Score2', 'CPICPhenotype', 'PAnnoPhenotype', 'Drug', 'Phenotypes', 'EvidenceLevel', 'LevelOverride', 'LevelModifier', 'Score', 'PMIDCount', 'EvidenceCount', 'Specialty', 'PhenotypeCategory']
ClinAnn = clin_ann_merge[table_cols].reset_index(drop = True)
ClinAnn.insert(0, 'ID', (pd.DataFrame(ClinAnn.index)+1)[0].to_list())

#columns_length(ClinAnn)
table_p = TSV["ClinAnn"]
ClinAnn.to_csv(table_p, sep="\t", index=0, quoting=3)

# Manual checking

# # Replace PAnno phenotype with CPIC phenotype
# ClinAnn = pd.read_csv(table_p, sep="\t")
# ClinAnn[ClinAnn.Score2.isnull()]
# response_score = {'Unknown': np.nan, 'Uncertain': np.nan, 'No': 0, 'Decreased': 0.5, 'Normal': 1, 'Increased': 2}

row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnn FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\r\\n' IGNORE 1 LINES (ID, CAID, Gene, Variant, Allele1, Allele2, Annotation1, Annotation2, Function1, Function2, @Score1, @Score2, CPICPhenotype, PAnnoPhenotype, Drug, Phenotypes, EvidenceLevel, LevelOverride, LevelModifier, Score, PMIDCount, EvidenceCount, Specialty, PhenotypeCategory) \
    SET Score1 = CASE WHEN @Score1 NOT IN ('None', '') THEN @Score1 END, \
      Score2 = CASE WHEN @Score2 NOT IN ('None', '') THEN @Score2 END;" % table_p)
row_count


#-------------PharmGKB.ClinAnnHasPhenotype-------------#
ClinAnnHasPhenotype = ClinAnn[["ID", "Phenotypes"]]
ClinAnnHasPhenotype = ClinAnnHasPhenotype.drop(['Phenotypes'], axis=1).join(ClinAnnHasPhenotype['Phenotypes'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Phenotype')).reset_index(drop = True)

# 注意：这里使用 - 来替代空值，大概会有两万多行-
ClinAnnHasPhenotype.Phenotype = ClinAnnHasPhenotype.Phenotype.replace(np.nan, "-")

# Separate access the PhenotypeID and then merge the two dataframes
phenotype_df = ClinAnnHasPhenotype[["Phenotype"]].drop_duplicates().reset_index(drop = True)
phenotype_df["PhenotypeID"] = ""
for index, row in phenotype_df.iterrows():
  row["PhenotypeID"] = pymysql_cursor("SELECT ID FROM Phenotype WHERE Name ='%s';" % row["Phenotype"])
  if type(row["PhenotypeID"]) is list:
    row["PhenotypeID"] = row["PhenotypeID"][0]  # 两次更新中都有'hypersexuality state'
  phenotype_df.iloc[index] = row

ClinAnnHasPhenotype = pd.merge(ClinAnnHasPhenotype, phenotype_df, on="Phenotype")
ClinAnnHasPhenotype = ClinAnnHasPhenotype[["ID", "PhenotypeID"]]
ClinAnnHasPhenotype = ClinAnnHasPhenotype[ClinAnnHasPhenotype["PhenotypeID"].isna() == False]

table_p = TSV["ClinAnnHasPhenotype"]
ClinAnnHasPhenotype.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnnHasPhenotype FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (ClinAnnID, PhenotypeID);" % table_p)
row_count


#-------------PharmGKB.ClinAnnHasVariant-------------#
ClinAnnHasVariant = ClinAnn[["ID", "Variant"]]
ClinAnnHasVariant = ClinAnnHasVariant.drop(['Variant'], axis=1).join(ClinAnnHasVariant['Variant'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('Variant')).reset_index(drop = True)

# Separate access the VariantID and then merge the two dataframes
Variant_df = ClinAnnHasVariant[["Variant"]].drop_duplicates().reset_index(drop = True)
Variant_df["VariantID"] = ""
for index, row in Variant_df.iterrows():
  row["VariantID"] = pymysql_cursor("SELECT ID FROM Variant WHERE Name ='%s';" % row["Variant"])
  if type(row["VariantID"]) is list:
    print(row["VariantID"])
    # row["VariantID"] = row["VariantID"][0]
  Variant_df.iloc[index] = row

ClinAnnHasVariant = pd.merge(ClinAnnHasVariant, Variant_df, on="Variant")
ClinAnnHasVariant = ClinAnnHasVariant[["ID", "VariantID"]]
ClinAnnHasVariant = ClinAnnHasVariant[ClinAnnHasVariant["VariantID"].isna() == False]

table_p = TSV["ClinAnnHasVariant"]
ClinAnnHasVariant.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnnHasVariant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (ClinAnnID, VariantID);" % table_p)
row_count


#-------------PharmGKB.ClinAnnHasDrug-------------#
ClinAnnHasDrug = ClinAnn[["ID", "Drug"]]

# Separate access the DrugID and then merge the two dataframes
Drug_df = ClinAnnHasDrug[["Drug"]].drop_duplicates().reset_index(drop = True)
Drug_df["DrugID"] = ""
for index, row in Drug_df.iterrows():
  row["DrugID"] = pymysql_cursor("SELECT ID FROM Drug WHERE Name ='%s';" % row["Drug"])
  if type(row["DrugID"]) is list:
    print(row["DrugID"])
    # row["DrugID"] = row["DrugID"][0]
  Drug_df.iloc[index] = row

ClinAnnHasDrug = pd.merge(ClinAnnHasDrug, Drug_df, on="Drug")
ClinAnnHasDrug = ClinAnnHasDrug[["ID", "DrugID"]]
ClinAnnHasDrug = ClinAnnHasDrug[ClinAnnHasDrug["DrugID"].isna() == False]

table_p = TSV["ClinAnnHasDrug"]
ClinAnnHasDrug.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnnHasDrug FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (ClinAnnID, DrugID);" % table_p)
row_count



#-------------PharmGKB.Guideline-------------#
#-------------PharmGKB.GuidelineDetail-------------#
import glob, json
from bs4 import BeautifulSoup
files = glob.glob("%s/*.json" % RAW['guideline_path'])
guideline = []
guideline_detail = []
for f in files:
  Source = re.sub(RAW['guideline_path'], '', f).split("_")[2]
  with open(f) as load_f:
    content = json.load(load_f)
    PAID = content["guideline"]["id"] # https://www.pharmgkb.org/guidelineAnnotation/PA166262321
    # Summary
    raw_sum = content["guideline"]["summaryMarkdown"]["html"]
    Summary = BeautifulSoup(raw_sum, 'html.parser').find('p').text.replace('\n', ' ')
    # PA166202781 is removed manually
    if re.findall('no recommendation|no dosing recommendation|no longer dosing recommendation|no gene-drug interaction|No action is needed for this gene-drug interaction|preclude making recommendations', Summary) == []:
      # Alternate info
      Alternate = 1 if content["guideline"]["alternateDrugAvailable"] else 0
      # Dosing info
      Dosing = 1 if content["guideline"]["dosingInformation"] else 0
      # Related Drugs and Genes
      drugs = content["guideline"]["relatedChemicals"]
      genes = content["guideline"]["relatedGenes"]
      for drug in drugs:
        Drug = drug['name']
        DrugID = pymysql_cursor("SELECT ID FROM Drug WHERE PAID ='%s';" % drug["id"])
        for gene in genes:
          Gene = gene['symbol']
          GeneID = pymysql_cursor("SELECT ID FROM Gene WHERE PAID ='%s';" % gene["id"])
          guideline.append([Source, PAID, Summary, Gene, Drug, GeneID, DrugID, Alternate, Dosing])
      
      # Recommendation
      html_str = content["guideline"]["textMarkdown"]["html"]
      soup = BeautifulSoup(html_str, 'html.parser')
      tables = soup.find_all('table')
      for table in tables:
        dat = BeautifulSoup(re.sub('\x01', '', re.sub('<br/>', '}{', re.sub('( )?<sup>[a-z, ]+</sup>', '', str(table)))), 'html.parser').text.strip()
        th = dat.split('\n\n\n\n\n')[0]
        if re.findall('recommendation', th, re.IGNORECASE):
          ths = th.split('\n')
          # print(PAID, ths)
          trs = []
          for tr in dat.split('\n\n\n\n\n')[1].split('\n\n'):
            if tr != '':
              trs.append(tr.strip().split('\n'))
          # Select phenotype and recommendations
          df = pd.DataFrame(trs, columns=ths)
          phenotype = df.iloc[:,0].to_list()
          recommendation = []
          if 'Recommendation Indications OTHER than liver transplantation' in df.columns:
            for j in range(len(phenotype)):
              recommendation.append('Liver transplantation: %s Other indications: %s' % (df.loc[j, 'Recommendation LIVER transplantation'], df.loc[j, 'Recommendation Indications OTHER than liver transplantation']))
          elif 'Implications for PEG-IFN alpha and RBV' in df.columns:
            for j in range(len(phenotype)):
              recommendation.append('Implications for PEG-IFN alpha and RBV: %s Implications for protease inhibitor combinations with PEG-IFN alpha and RBV therapy: %s' % (df.loc[j, 'Implications for PEG-IFN alpha and RBV'], df.loc[j, 'Implications for protease inhibitor combinations with PEG-IFN alpha and RBV therapy']))
          elif 'Therapeutic recommendation - AIs are viable treatment option' in df.columns:
            for j in range(len(phenotype)):
              recommendation.append('%s or %s when AIs are contraindicated' % (df.loc[j, 'Therapeutic recommendation - AIs are viable treatment option'].replace('AI', 'aromatase inhibitor (AI)'), df.loc[j, 'Therapeutic recommendation - AIs are contraindicated']))
          else:
            sub_df = df.iloc[:, df.columns.str.match("^Recommendation(s)?$|^Therapeutic Dose Recommendation(s)?$|^Dosing recommendations.*$|^Recommendations for.*$|^Therapeutic Recommendation(s)?$|^RECOMMENDATION$", case=False)]
            recommendation = sub_df.iloc[:,0].to_list()
          if len(recommendation) != len(phenotype):
            print(df)
          for k in range(len(phenotype)):
            if len(genes) == 1 and len(drugs) == 1:
              guideline_detail.append([PAID, phenotype[k], recommendation[k], genes[0]['symbol'], drugs[0]['name']])
            else:
              guideline_detail.append([PAID, phenotype[k], recommendation[k], 'X', 'X'])

Guideline = pd.DataFrame(guideline, columns=["Source", "PAID", "Summary", "Gene", "Drug", "GeneID", "DrugID", "Alternate", "Dosing"])
Guideline.to_csv(TSV["Guideline"], sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Guideline FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Source, PAID, Summary, Gene, Drug, GeneID, DrugID, Alternate, Dosing);" % TSV["Guideline"]); row_count

GuidelineDetail = pd.DataFrame(guideline_detail, columns=["PAID", "Phenotype", "Recommendation", "Gene", "Drug"])
GuidelineDetail.to_csv(TSV["GuidelineDetail"], sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE GuidelineDetail FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (PAID, Phenotype, Genotype, Recommendation, Avoid, Gene, Drug);" % TSV["GuidelineDetail"]); row_count


#-------------PharmGKB.GuidelineMerge-------------#
### Mapping the Variant/Allele1,2
Guideline = pd.read_csv(TSV["Guideline"], sep="\t")
GuidelineDetail = pd.read_csv(TSV["GuidelineDetail"], sep="\t")
GuidelineMerge = Guideline.merge(GuidelineDetail, how='left', on='PAID')
# Check duplicated
GuidelineMerge[GuidelineMerge.duplicated()]
# Check drug and gene
tmp_index = []
for index, row in GuidelineMerge.iterrows():
  if type(row.Drug_y) is str:
    if row.PAID == 'PA166254881':
      row.Drug_y = row.Drug_x
    if row.Gene_x != row.Gene_y or row.Drug_x != row.Drug_y:
      tmp_index.append(index)
  else:
    row.Drug_y = row.Drug_x
    row.Gene_y = row.Gene_x
  GuidelineMerge.iloc[index] = row
# Manual check
GuidelineMerge.iloc[tmp_index,].to_csv('~/Downloads/GuidelineMerge.txt', sep='\t')
GuidelineMerge = GuidelineMerge.drop(columns=['Gene_x', 'Drug_x']).rename(columns={'Gene_y':'Gene', 'Drug_y':'Drug'})
GuidelineMerge.to_csv('./data/pharmgkb/intermediate/GuidelineMerge.txt', sep="\t", index=0, quoting=3)


GuidelineMerge = pd.read_csv('./data/pharmgkb/intermediate/GuidelineMerge.txt', sep="\t")
# Split Drugs
GuidelineMerge = GuidelineMerge.drop(['Drug'], axis=1).join(GuidelineMerge['Drug'].str.split('; ', expand=True).stack().reset_index(level=1, drop=True).rename('Drug')).drop_duplicates().reset_index(drop = True)
for index, row in GuidelineMerge.iterrows():
  # test = pymysql_cursor("SELECT ID FROM Guideline WHERE Gene ='%s' AND PAID='%s';" % (row.Gene, row.PAID))
  # if test is None:
  #   print(index)
  row.GeneID = pymysql_cursor("SELECT ID FROM Gene WHERE Symbol ='%s';" % row.Gene)
  row.DrugID = pymysql_cursor("SELECT ID FROM Drug WHERE Name ='%s';" % row.Drug)
  if row.DrugID is None:
    print(row)
  GuidelineMerge.iloc[index] = row


### Inserting the empty manually
GuidelineMerge = GuidelineMerge.drop_duplicates()
GuidelineMerge.insert(0, 'ID', (pd.DataFrame(GuidelineMerge.index)+1)[0].to_list())
GuidelineMerge = GuidelineMerge[['ID', 'Source', 'PAID', 'Summary', 'Phenotype', 'Genotype', 'Recommendation', 'Avoid', 'Alternate', 'Dosing', 'Gene', 'Drug', 'GeneID', 'DrugID']]
GuidelineMerge.to_csv(TSV["GuidelineMerge"], sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE GuidelineMerge FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (ID, Source, PAID, Summary, Phenotype, Genotype, Recommendation, Avoid, Alternate, Dosing, Gene, Drug, GeneID, DrugID);" % TSV["GuidelineMerge"]); row_count



#-------------PharmGKB.GuidelineRule-------------#
GuidelineMerge = pd.read_csv(TSV["GuidelineMerge"], sep="\t")
GuidelineMerge = GuidelineMerge[GuidelineMerge.ID >= 632]

# GuidelineMerge = GuidelineMerge.drop(columns=['ID']).drop_duplicates().reset_index(drop = True)
# GuidelineMerge.insert(0, 'ID', (pd.DataFrame(GuidelineMerge.index)+1)[0].to_list())
# GuidelineMerge = GuidelineMerge[['ID', 'Source', 'PAID', 'Summary', 'Phenotype', 'Genotype', 'Recommendation', 'Avoid', 'Alternate', 'Dosing', 'Gene', 'Drug', 'GeneID', 'DrugID']]
# GuidelineMerge.to_csv(TSV["GuidelineMerge"], sep="\t", index=0, quoting=3)

GuidelineMerge = GuidelineMerge.drop(['Genotype'], axis=1).join(GuidelineMerge['Genotype'].str.split('; ', expand=True).stack().reset_index(level=1, drop=True).rename('Genotype')).drop_duplicates().reset_index(drop = True)
GuidelineMerge = GuidelineMerge.drop(['Phenotype'], axis=1).join(GuidelineMerge['Phenotype'].str.split('; | \+ ', expand=True).stack().reset_index(level=1, drop=True).rename('Phenotype')).drop_duplicates().reset_index(drop = True)



guideline_rule = []
for index, row in GuidelineMerge.iterrows():
  print("%s/%s" % (index, GuidelineMerge.shape[0]))
  if row.Drug in ['daunorubicin', 'doxorubicin']:
    Drug = 'anthracyclines and related substances'
  else:
    Drug = row.Drug
  gt = row.Genotype
  if type(gt) is str:
    # SNP/Indel
    if gt.startswith('rs'):
      Variant = gt.split(' ')[0]
      if len(gt.split(' ')) == 2:
        alleles = gt.split(' ')[1]
        if '/' in alleles:
          Allele1 = alleles.split('/')[0]; Allele2 = alleles.split('/')[1]
        else:
          Allele1 = alleles[0]; Allele2 = alleles[1]
          ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Variant = "%s" AND Allele1 = "%s" AND Allele2 = "%s" AND Drug = "%s";' % (row.Gene, Variant, Allele1, Allele2, Drug))
          if ClinAnnID is None:
            ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Variant = "%s" AND Allele1 = "%s" AND Allele2 = "%s" AND Drug = "%s";' % (row.Gene, Variant, Allele2, Allele1, Drug))
            if ClinAnnID is None:
              continue # print(row.Gene, Variant, alleles, Drug) # UGT1A6 rs17863783 TT anthracyclines and related substances is Level 4, not included.
      else:
        ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Variant = "%s" AND Drug = "%s";' % (row.Gene, Variant, Drug))
        if ClinAnnID is None:
          continue # print(row.Gene, Variant, alleles, Drug) # VKORC1 rs9923231 AG fluindione
      # Collect the results
      if type(ClinAnnID) is list:
        ClinAnnID_str = '; '.join([str(s) for s in ClinAnnID])
      else:
        ClinAnnID_str = ClinAnnID
      guideline_rule.append([row.Gene, Variant, Allele1, Allele2, row.Phenotype, ClinAnnID_str, row.ID])
    # Haplotypes
    elif gt.startswith('AS'):
      Variant = ''
      res = cursor.execute('SELECT Allele1, Allele2 FROM DiplotypePhenotype WHERE Gene = "%s" AND ActivityScore = "%s";' % (row.Gene, format(float(gt.split(' ')[1]), ".1f")))
      res = cursor.fetchall()
      if res == ():
        print(row)
      else:
        for alleles in res:
          Allele1 = alleles[0]; Allele2 = alleles[1]
          ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s" AND Drug = "%s";' % (row.Gene, Allele1, Allele2, Drug))
          if ClinAnnID is None:
            ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s" AND Drug = "%s";' % (row.Gene, Allele2, Allele1, Drug))
            if ClinAnnID is None:
              ClinAnnID = ''
              # print(row.Gene, alleles, Drug) # DPYD and CYP2D6
              # pymysql_cursor('SELECT Variant FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (row.Gene, Allele1))
          # Collect the results
          if type(ClinAnnID) is list:
            ClinAnnID_str = '; '.join([str(s) for s in ClinAnnID])
          else:
            ClinAnnID_str = ClinAnnID
          guideline_rule.append([row.Gene, Variant, Allele1, Allele2, row.Phenotype, ClinAnnID_str, row.ID])
    elif gt.startswith('*'):
      Variant = ''
      if '/' in gt:
        alleles = gt.split('/')
        Allele1 = alleles[0]; Allele2 = alleles[1]
      else:
        Allele1 = gt; Allele2 = ''
      # Collect the results
      guideline_rule.append([row.Gene, Variant, Allele1, Allele2, row.Phenotype, '', row.ID])
    else: # CFTR
      Variant = pymysql_cursor('SELECT Variant FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (row.Gene, gt))
      Allele1 = pymysql_cursor('SELECT AlleleManual FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (row.Gene, gt))
      Allele2 = ''
      # Collect the results
      guideline_rule.append([row.Gene, Variant, Allele1, Allele2, row.Phenotype, '', row.ID])
  else:
    phe = row.Phenotype
    if phe.startswith('HLA'):
      Variant = phe
      Allele1 = phe.split(' ')[0].replace(row.Gene, '')
      Allele2 = ''
      # Collect the results
      guideline_rule.append([row.Gene, Variant, Allele1, Allele2, row.Phenotype, '', row.ID])
    else:
      Variant = ''
      phe = phe.replace('CYP2C9 ', '')
      res = cursor.execute('SELECT Allele1, Allele2 FROM DiplotypePhenotype WHERE Gene = "%s" AND Phenotype = "%s";' % (row.Gene, phe))
      res = cursor.fetchall()
      if phe == 'Intermediate Metabolizer - *10':
        res = cursor.execute('SELECT Allele1, Allele2 FROM DiplotypePhenotype WHERE Gene = "%s" AND Phenotype = "%s" AND Allele1 != "*10" AND Allele2 != "*10" AND Allele1 != "*10x2" AND Allele2 != "*10x2" AND ActivityScore = 1;' % (row.Gene, phe))
        res = cursor.fetchall()
      if phe == 'Intermediate Metabolizer' and row.PAID == 'PA166181885':
        res = cursor.execute('SELECT Allele1, Allele2 FROM DiplotypePhenotype WHERE Gene = "%s" AND Phenotype = "%s" AND ID NOT IN (SELECT ID FROM DiplotypePhenotype WHERE Gene = "%s" AND Phenotype = "%s" AND Allele1 != "*10" AND Allele2 != "*10" AND Allele1 != "*10x2" AND Allele2 != "*10x2" AND ActivityScore = 1);' % (row.Gene, phe, row.Gene, phe))
        res = cursor.fetchall()
      for alleles in res:
        Allele1 = alleles[0]; Allele2 = alleles[1]
        ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s" AND Drug = "%s";' % (row.Gene, Allele1, Allele2, Drug))
        if ClinAnnID is None:
          ClinAnnID = pymysql_cursor('SELECT ID FROM ClinAnn WHERE Gene = "%s" AND Allele1 = "%s" AND Allele2 = "%s" AND Drug = "%s";' % (row.Gene, Allele2, Allele1, Drug))
          if ClinAnnID is None:
            ClinAnnID = ''
            # print(row.Gene, alleles, Drug) # DPYD and CYP2D6
            # pymysql_cursor('SELECT Variant FROM AlleleFunctionality WHERE Gene = "%s" AND Allele = "%s";' % (row.Gene, Allele1))
        # Collect the results
        if type(ClinAnnID) is list:
          ClinAnnID_str = '; '.join([str(s) for s in ClinAnnID])
        else:
          ClinAnnID_str = ClinAnnID
        guideline_rule.append([row.Gene, Variant, Allele1, Allele2, row.Phenotype, ClinAnnID_str, row.ID])

GuidelineRule = pd.DataFrame(guideline_rule, columns=['Gene', 'Variant', 'Allele1', 'Allele2', 'Phenotype', 'ClinAnnID', 'GuidelineID']).drop_duplicates()
columns_length(GuidelineRule)
GuidelineRule.to_csv(TSV["GuidelineRule"], sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE GuidelineRule FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Gene, Variant, Allele1, Allele2, Phenotype, ClinAnnID, GuidelineID);" % TSV["GuidelineRule"]); row_count

# Manual update
# update `PAnno`.`GuidelineRule` set Phenotype='HLA-B*57:01 positive' WHERE Phenotype='HLA-B*5701';
# update `PAnno`.`GuidelineRule` set Variant='HLA-B*57:01 positive' WHERE Variant='HLA-B*5701';
# update `PAnno`.`GuidelineRule` set Allele1='*57:01' WHERE Allele1='*5701';
# update `PAnno`.`GuidelineRule` set Allele1='*15:02' WHERE Allele1='HLA-B*15:02';
# update `PAnno`.`GuidelineRule` set Allele1='*31:01' WHERE Allele1='HLA-A*31:01';
# update `PAnno`.`GuidelineRule` set Phenotype='Poor Metabolizer' WHERE Phenotype='CYP2C9 Poor Metabolizer';
# update `PAnno`.`GuidelineRule` set Phenotype='Intermediate Metabolizer' WHERE Phenotype='CYP2C9 Intermediate Metabolizer';
# update `PAnno`.`GuidelineRule` set Phenotype='Normal Metabolizer' WHERE Phenotype='CYP2C9 Normal Metabolizer';
# update `PAnno`.`GuidelineRule` set Phenotype='Decreased Function' WHERE Gene='SLCO1B1' AND Allele1='*1' AND Allele2='*2';
# update `PAnno`.`GuidelineRule` set Phenotype='Intermediate Metabolizer' WHERE Gene='CYP2B6' AND Phenotype='*5/*6 or *5/*18';



# ### Dump the final version and check manually
# guide = cursor.execute("SELECT * FROM GuidelineMerge")
# guide = cursor.fetchall()
# guide_df = pd.DataFrame(guide, columns=['ID', 'Source', 'PAID', 'Summary', 'Phenotype', 'Genotype', 'Recommendation', 'Avoid', 'Alternate', 'Dosing', 'Gene', 'Drug', 'GeneID', 'DrugID'])
# guide_df.Summary = guide_df.Summary.apply(lambda x : x.replace('\n', '\\n'))
# guide_df.Recommendation = guide_df.Recommendation.apply(lambda x: x.replace('\n \n', '\n')).apply(lambda x : x.replace('\n', '\\n'))
# guide_df.to_csv(TSV["GuidelineMerge"], sep="\t", index=0)#, quoting=3)


# ### Dump the final version and check manually
# rule = cursor.execute("SELECT * FROM GuidelineRule")
# rule = cursor.fetchall()
# rule_df = pd.DataFrame(rule, columns=["ID", "Gene", "Variant", "Allele1", "Allele2", "Phenotype", "ClinAnnID", "GuidelineID"])
# rule_df = rule_df.drop(columns = ["ID"])
# rule_df.to_csv('./data/panno/tables/PharmGKB.GuidelineRule.txt', sep="\t", index=0)#, quoting=3)
