#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Process PharmGKB data.
"""

# Import modules
import os, json, html, toml, re
import pandas as pd
import numpy as np
from utils.pgi_functions import connect_database, pymysql_cursor, load_logging_cfg
from modules.pharmgkb.parse_anno import process_uncomplete

cfile = './conf/pharmgkb.toml'
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
#-------------PharmGKB.Gene-------------#
# Load data
genes = pd.read_csv(RAW["genes"], sep="\t", dtype="object")
#ann_genes = genes[(genes["Is VIP"] == "Yes") | (genes["Has Variant Annotation"] == "Yes") | (genes["Has CPIC Dosing Guideline"] == "Yes")]
Gene = genes[['Symbol', 'Name', 'PharmGKB Accession Id', 'NCBI Gene ID', 'HGNC ID', 'Ensembl Id', 'Chromosome', 'Chromosomal Start - GRCh37', 'Chromosomal Stop - GRCh37', 'Chromosomal Start - GRCh38', 'Chromosomal Stop - GRCh38']].reset_index(drop = True)

table_p = TSV["Gene"]
Gene.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Gene FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Symbol, Name, PAID, EntrezID, HGNCID, EnsemblID, Chr, GRCh37Start, GRCh37End, GRCh38Start, GRCh38End);" % table_p)
row_count

# Insert the '-'
pymysql_cursor('INSERT INTO Gene(Symbol) VALUES("-");')


#-------------PharmGKB.Variant-------------#
# Load data
variants = pd.read_csv(RAW["variants"], sep="\t", dtype="object")
Variant = variants[['Variant Name', 'Location', 'Variant ID']]
Variant["Type"] = "Variant"

# Add Haplotype into PharmGKB.Variant table as a specific variant
# Load data
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
VariantOrHaplotype.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Variant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name, NCID, PAID, Type);" % table_p)
row_count


#-------------PharmGKB.VariantSynonyms-------------#
# Load data
variants = pd.read_csv(RAW["variants"], sep="\t")
VariantSynonyms = variants[['Variant ID', 'Synonyms']].reset_index(drop = True)
for index, row in VariantSynonyms.iterrows():
  VariantSynonyms.loc[index, "Variant ID"] = pymysql_cursor("SELECT ID FROM Variant Where PAID = '%s';" % row["Variant ID"])

VariantSynonyms = VariantSynonyms.drop(['Synonyms'], axis=1).join(VariantSynonyms['Synonyms'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('Synonym'))
columns_length(VariantSynonyms)

table_p = TSV["VariantSynonyms"]
VariantSynonyms.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE VariantSynonyms FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (VariantID, Synonym);" % table_p)
row_count


#-------------PharmGKB.GeneHasVariant-------------#
# Load data
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
GeneHasVariantOrHaplotype.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE GeneHasVariant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (VariantID, GeneID);" % table_p)
row_count


#-------------PharmGKB.Drug-------------#
# Load data
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
Drug.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Drug FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ID, Name, PAID, GenericName, TradeName, CrossRef, Type);" % table_p)
row_count


#-------------PharmGKB.DrugSynonyms-------------#
# Load data
DrugSynonyms = DrugSynonyms.dropna().reset_index(drop = True)
DrugSynonyms = DrugSynonyms.drop(['Synonyms'], axis=1).join(DrugSynonyms['Synonyms'].str.split(',"', expand=True).stack().reset_index(level=1, drop=True).rename('Synonym').str.strip('"'))
DrugSynonyms = DrugSynonyms.drop_duplicates().reset_index(drop = True)

DrugSynonyms = DrugSynonyms[DrugSynonyms["Synonym"] != ""]
columns_length(DrugSynonyms)

table_p = TSV["DrugSynonyms"]
DrugSynonyms.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE DrugSynonyms FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (DrugID, Synonym);" % table_p)
row_count


#-------------PharmGKB.Phenotype-------------#
# Load data
phenotypes = pd.read_csv(RAW["phenotypes"], sep="\t")
Phenotype = phenotypes[["PharmGKB Accession Id", "Name", "Alternate Names"]].drop_duplicates()
Phenotype.insert(0, 'ID', (pd.DataFrame(Phenotype.index)+1)[0].to_list())

table_p = TSV["Phenotype"]
Phenotype[["ID", "PharmGKB Accession Id", "Name"]].to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Phenotype FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ID, PAID, Name);" % table_p)
row_count

pymysql_cursor('INSERT INTO Phenotype(Name, PAID, Type) VALUES("-", "-", "-");')


#-------------PharmGKB.PhenotypeSynonyms-------------#
# Load data
Phenotype[["Name", "Alternate Names"]] = Phenotype[["Name", "Alternate Names"]].replace(np.nan, "")
Phenotype["Alternate Names"] = Phenotype["Name"].str.cat(Phenotype["Alternate Names"], sep = ',"')

PhenotypeSynonyms = Phenotype[["ID", "Alternate Names"]].drop_duplicates().rename(columns = {'ID': 'PhenotypeID'}).reset_index(drop = True)
PhenotypeSynonyms = PhenotypeSynonyms.drop(['Alternate Names'], axis=1).join(PhenotypeSynonyms['Alternate Names'].str.split(',"', expand=True).stack().reset_index(level=1, drop=True).rename('Synonym').str.strip('"'))
PhenotypeSynonyms = PhenotypeSynonyms.drop_duplicates()
PhenotypeSynonyms = PhenotypeSynonyms[PhenotypeSynonyms["Synonym"] != ""]
columns_length(PhenotypeSynonyms)

table_p = TSV["PhenotypeSynonyms"]
PhenotypeSynonyms.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PhenotypeSynonyms FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (PhenotypeID, Synonym);" % table_p)
row_count

#### Determine if this phenotype is cancer
cancer_reg = "|".join(["cancer", "cancers", "tumors", "tumor", "tumours", "tumour", "carcinomas", "carcinoma", "neoplasms", "neoplasm", "epitheliomas", "epithelioma", "malignancies", "malignancy", "adenocarcinoma", "cystadenoma", "cystadenocarcinoma", "sarcoma"])
for index, row in PhenotypeSynonyms.iterrows():
  if re.findall(cancer_reg, row.Synonym.lower()):
    pymysql_cursor("UPDATE PharmGKB.Phenotype SET Type = 'Cancer' WHERE ID = '%s';" % row.PhenotypeID)
  else:
    pymysql_cursor("UPDATE PharmGKB.Phenotype SET Type = 'Not cancer' WHERE ID = '%s';" % row.PhenotypeID)



#########################################
## Annotation Related Dictionary
#########################################
#-------------PharmGKB.EvidenceTypeDic-------------#
# Load data
clinical_ann_evidence = pd.read_csv(RAW["clinical_ann_evidence"], sep="\t")
EvidenceTypeDic = pd.DataFrame({"Name": list(set(clinical_ann_evidence["Evidence Type"].to_list()))})

table_p = TSV["EvidenceTypeDic"]
EvidenceTypeDic.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE EvidenceTypeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % table_p)
row_count


#-------------PharmGKB.PhenotypeCategoryDic-------------#
clinical_annotations = pd.read_csv(RAW["clinical_annotations"], sep="\t")
pheno_category = list(set(clinical_annotations["Phenotype Category"].str.split(";", expand=True).stack().to_list())) + ["PD"]
PhenotypeCategoryDic = pd.DataFrame({"Name": pheno_category})

table_p = TSV["PhenotypeCategoryDic"]
PhenotypeCategoryDic.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PhenotypeCategoryDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % table_p)
row_count




#########################################
## Variant Annotation
#########################################
#-------------PharmGKB.VarAnn-------------#
# Load data
# Replace CYP2C19􏰀17 into CYP2C19*17
# 􏰃1.569 to <2.075 mg/kg/ day
# 􏰃1.54 to <2.05 mg/kg/day
var_drug_ann = pd.read_csv(RAW["var_drug_ann"], sep="\t")
var_fa_ann = pd.read_csv(RAW["var_fa_ann"], sep="\t")
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
VarAnn.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor('LOAD DATA LOCAL INFILE "%s" INTO TABLE VarAnn FIELDS TERMINATED BY "\\t" LINES TERMINATED BY "\\n" (AnnotationID, VariantOrHaplotype, Gene, Drug, PMID, Significance, Note, Sentence, Allele, SpecialtyPopulation, TypeID, PhenotypeCategoryID);' % table_p)
row_count




#########################################
## Clinical Annotation
#########################################
#-------------PharmGKB.ClinAnnEvidence-------------#
# Load data
clinical_ann_evidence = pd.read_csv(RAW["clinical_ann_evidence"], sep="\t", dtype="object")

for index, row in clinical_ann_evidence.iterrows():
  row["Evidence Type"] = pymysql_cursor("SELECT ID FROM EvidenceTypeDic WHERE Name ='%s';" % row["Evidence Type"])
  clinical_ann_evidence.iloc[index] = row

clinical_ann_evidence["PMID"] = clinical_ann_evidence["PMID"].replace(np.nan, '')
table_p = TSV["ClinAnnEvidence"]
clinical_ann_evidence.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnnEvidence FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (CAID, EvidenceID, EvidenceTypeID, URL, PMID, Summary, Score);" % table_p)
row_count


#-------------PharmGKB.ClinAnn-------------#
# Load data
clinical_annotations = pd.read_csv(RAW["clinical_annotations"], sep="\t")
clinical_ann_alleles = pd.read_csv(RAW["clinical_ann_alleles"], sep="\t")
# Merge annotations and alleles
clinical_ann = pd.merge(clinical_annotations, clinical_ann_alleles, on="Clinical Annotation ID")

# Replace the nan into -
nan_index = clinical_ann[clinical_ann['Gene'].isna()].index.to_list()
clinical_ann.loc[nan_index, 'Gene'] = "-"
nan_index = clinical_ann[clinical_ann['Phenotype Category'].isna()].index.to_list()
clinical_ann.loc[nan_index,'Phenotype Category'] = "Other"

# Splite the variants and phenotype categories
clinical_ann = clinical_ann.drop(['Variant/Haplotypes'], axis=1).join(clinical_ann['Variant/Haplotypes'].str.split(', ', expand=True).stack().reset_index(level=1, drop=True).rename('Variant/Haplotype')).reset_index(drop = True)
clinical_ann = clinical_ann.drop(['Phenotype Category'], axis=1).join(clinical_ann['Phenotype Category'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Phenotype Category')).reset_index(drop = True)

# Remove some duplications happened when mergeing. For example, clinical_ann[clinical_ann['Clinical Annotation ID'] == 1184746746]
retained_index = []
for index, row in clinical_ann.iterrows():
  gene = row['Gene']
  alleles = row['Genotype/Allele']
  # Cleaning data
  if type(alleles) == str:
    alleles = re.sub('A-968C_376G', 'A- 968C_376G', re.sub('A-680T_376G', 'A- 680T_376G', re.sub('A-202A_376G', 'A- 202A_376G', alleles)))
    clinical_ann.loc[index, 'Genotype/Allele'] = alleles
  haplotype = row['Variant/Haplotype'].replace(gene, '').strip()
  # Single locus
  if haplotype.startswith('rs') or type(alleles) != str:
    retained_index.append(index)
  # Others
  elif haplotype in alleles:
    retained_index.append(index)

clinical_ann.shape; len(retained_index)
clinical_ann = clinical_ann.iloc[retained_index,].reset_index(drop = True)
clinical_ann.shape

####### BE ATTENTION #######
# Split the clinical_ann into two parts
clinical_ann["Comment"] = ""
# # Only extremely clear functions are retained.
# # Only retain 'Normal', 'Increased', 'No', 'Decreased', 'Uncertain', 'Unknown'
# nan_rows = clinical_ann[clinical_ann["Allele Function"].isna()]
# unnan_rows = clinical_ann[clinical_ann["Allele Function"].isna() == False]
# unnan_rows["Allele Function"] = unnan_rows["Allele Function"].str.split(' ', expand=True)[0]
# unnan_rows["Allele Function"].drop_duplicates()
# uncom = unnan_rows[unnan_rows["Allele Function"].str.contains('Normal|Increased|No|Decreased|Uncertain|Unknown', regex=True) == False]

# # Uncomplete and Complete
# uncomplete = pd.concat([nan_rows, uncom]).reset_index(drop = True)
# complete = unnan_rows[unnan_rows["Allele Function"].str.contains('Normal|Increased|No|Decreased|Uncertain|Unknown', regex=True)].reset_index(drop = True)
# clinical_ann.shape; uncomplete.shape; complete.shape
# uncomplete.to_csv("./data/pharmgkb/intermediate/uncomplete.txt", sep="\t", index=0)
# complete.to_csv("./data/pharmgkb/intermediate/complete.txt", sep="\t", index=0)

## 上面的部分都被注释掉了，因为网站自带的注释其实和最终的药效还是有一定的gap
## 例如: https://www.pharmgkb.org/clinicalAnnotation/1184746746同时包括了剂量与毒性
## TPMT*15被注释为No function，但是在剂量方面是decreased状态，在毒性方面是increased状态
need_manual, rule_res = process_uncomplete(clinical_ann)
need_manual.shape; rule_res.shape; clinical_ann.shape
need_manual.to_csv("./data/pharmgkb/intermediate/need_manual.txt", sep="\t", index=0)

##### Manual process
### reset_index 超级超级重要！！！！
manualed = pd.read_csv("./data/pharmgkb/intermediate/need_manual.txt", sep="\t")
all_clinical_ann = pd.concat([rule_res, manualed]).reset_index(drop = True)
all_clinical_ann.shape

############################
# Splite the drugs and genes
all_clinical_ann = all_clinical_ann.drop(['Drug(s)'], axis=1).join(all_clinical_ann['Drug(s)'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Drug')).reset_index(drop = True)
all_clinical_ann = all_clinical_ann.drop(['Gene'], axis=1).join(all_clinical_ann['Gene'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Gene')).reset_index(drop = True)
all_clinical_ann.shape

# DicID
all_clinical_ann = pd.concat([all_clinical_ann, pd.DataFrame(columns=["PhenotypeCategoryID", "GeneID", "VariantID", "DrugID"])], axis=1)
for index, row in all_clinical_ann.iterrows():
  # print("%s/%s" % (index, all_clinical_ann.shape[0]))
  row["PhenotypeCategoryID"] = pymysql_cursor("SELECT ID FROM PhenotypeCategoryDic WHERE Name ='%s';" % row["Phenotype Category"])
  row["GeneID"] = pymysql_cursor("SELECT ID FROM Gene WHERE Symbol ='%s';" % row["Gene"])
  row["VariantID"] = pymysql_cursor("SELECT ID FROM Variant WHERE Name ='%s';" % row["Variant/Haplotype"])
  row["DrugID"] = pymysql_cursor("SELECT ID FROM Drug WHERE Name ='%s';" % row["Drug"])
  if row["PhenotypeCategoryID"] is None or row["GeneID"] is None or row["VariantID"] is None or row["DrugID"] is None:
    print(index)
  all_clinical_ann.iloc[index] = row


# Insert ID
all_clinical_ann['Clinical Annotation ID'] = all_clinical_ann['Clinical Annotation ID'].astype('int')
all_clinical_ann = all_clinical_ann.reset_index(drop=True)
all_clinical_ann.insert(0, 'ID', (pd.DataFrame(all_clinical_ann.index)+1)[0].to_list())
ClinAnn = all_clinical_ann[['ID', 'Clinical Annotation ID', 'Gene', 'Variant/Haplotype', 'Drug', 'Phenotype(s)', 'Level of Evidence', 'Level Override', 'Level Modifiers', 'Score', 'PhenotypeCategoryID', 'Genotype/Allele', 'Annotation Text', 'Allele Function', 'PMID Count', 'Evidence Count', 'URL', 'Specialty Population', 'GeneID', 'VariantID', 'DrugID', 'Comment']]

#columns_length(ClinAnn)
table_p = TSV["ClinAnn"]
ClinAnn.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnn FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ID, CAID, Gene, VariantOrHaplotype, Drug, Phenotypes, EvidenceLevel, LevelOverride, LevelModifier, Score, PhenotypeCategoryID, GenotypeOrAllele, Annotation, Function, PMIDCount, EvidenceCount, URL, SpecialtyPopulation, GeneID, VariantID, DrugID, Comment);" % table_p)
row_count


#-------------PharmGKB.ClinAnnHasPhenotype-------------#
ClinAnnHasPhenotype = ClinAnn[["ID", "Phenotype(s)"]]
ClinAnnHasPhenotype = ClinAnnHasPhenotype.drop(['Phenotype(s)'], axis=1).join(ClinAnnHasPhenotype['Phenotype(s)'].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('Phenotype')).reset_index(drop = True)

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
ClinAnnHasPhenotype.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinAnnHasPhenotype FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, PhenotypeID);" % table_p)
row_count




#-------------PharmGKB.ClinDosingGuideline-------------#
files = os.listdir(RAW['guideline_path'])
guidelines = []
for f in files:
  if ".json" in f:
    Source = f.split("_")[2]
    with open("%s/%s" % (RAW['guideline_path'], f)) as load_f:
      content = json.load(load_f)
      URL = content["guideline"]["@id"]
      raw_sum = html.unescape(content["guideline"]["summaryMarkdown"]["html"]).strip()
      Annotation = re.sub('<em>|</em>|<p>|</p>|<strong>|</strong>|<a href=".*">|</a>|<a rel=".*">', '', raw_sum)
      Annotation = re.sub('\n', '', Annotation)
      # Remove '"", etc.
      Annotation = Annotation.strip('"').replace('""', "'").replace('"', "'")
      # Related Drugs and Genes
      drugs = content["guideline"]["relatedChemicals"]
      genes = content["guideline"]["relatedGenes"]
      for drug in drugs:
        RelatedDrugID = pymysql_cursor("SELECT ID FROM Drug WHERE PAID ='%s';" % drug["id"])
        for gene in genes:
          RelatedGeneID = pymysql_cursor("SELECT ID FROM Gene WHERE PAID ='%s';" % gene["id"])
          guidelines.append([Source, Annotation, RelatedGeneID, RelatedDrugID, URL])

ClinDosingGuideline = pd.DataFrame(guidelines, columns=["Source", "Annotation", "RelatedGeneID", "RelatedDrugID", "URL"])

table_p = TSV["ClinDosingGuideline"]
ClinDosingGuideline.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE ClinDosingGuideline FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Source, Annotation, RelatedGeneID, RelatedDrugID, URL);" % table_p)
row_count


"""
注意，这里不强行往Annotation中增加Dosing Guidelines了，而是改为：在向OT中插入PharmGKB的数据时，查询PharmGKB.ClinDosingGuidelines

## Add the guideline into ClinAnn
ClinDosingGuideline['ClinAnnIDs'] = ''
for index, row in ClinDosingGuideline.iterrows():
  # Get ClinAnnID
  ClinAnnIDs = ClinAnn[(ClinAnn.GeneID == row.RelatedGeneID) & (ClinAnn.DrugID == row.RelatedDrugID)].ID.to_list()
  row.ClinAnnIDs = ",".join([str(i) for i in ClinAnnIDs])
  ClinDosingGuideline.iloc[index] = row

map_df = []
for ele in set(ClinDosingGuideline.ClinAnnIDs):
  if len(ele) > 0:
    df = ClinDosingGuideline[ClinDosingGuideline.ClinAnnIDs == ele]
    merge = []
    for index, row in df.iterrows():
      source_ann = "%s: %s" % (row.Source, row.Annotation)
      merge.append(source_ann)
    map_df.append([ele, "; ".join(merge)])


for ele in map_df:
  ClinAnnIDs = ele[0].split(",")
  for caid in ClinAnnIDs:
    index = int(caid)-1
    rawann = ClinAnn.loc[index, "Annotation Text"]
    dosing_rawann = "Dosing Guideline: %s | Genotype Function: %s" % (ele[1], rawann)
    ClinAnn.loc[index, "Annotation Text"] = dosing_rawann
"""