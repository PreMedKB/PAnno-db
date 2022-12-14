#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Map PharmGKB primary data into MetaData.
"""

# Import modules
import toml
import re
import pandas as pd
from utils.pgi_functions import connect_database, pymysql_cursor, load_logging_cfg, map_meta
from utils.trim import drug_trim, disease_trim

cfile = './conf/pharmgkb.toml'
cf = toml.load(cfile, _dict=dict)
DEFAULT = cf['DEFAULT']
RAW = cf['RAW']
TSV = cf['TSV']

global db, cursor
db, cursor = connect_database(DEFAULT)
logger = load_logging_cfg(DEFAULT)


# Noted: The primary data of gene in PharmGKB exist some errors, such as the expired gene symbol and ID.
#        Therefore, we only map the genes which have clinical annotations into PreMedKB MetaData.

##### Gene
gene2meta  = []
Gene = cursor.execute('SELECT GeneID, Gene From ClinAnn;')
Gene = cursor.fetchall()
Gene = pd.DataFrame(Gene, columns = ["GeneID", "Symbol"]).drop_duplicates()

for index, row in Gene.iterrows():
  mg_id = pymysql_cursor('SELECT ID FROM Meta.Gene WHERE Symbol = "%s";' % row.Symbol)
  if mg_id is None:
    mg_id = pymysql_cursor('SELECT GeneID FROM Meta.GeneAlias WHERE Alias = "%s";' % row.Symbol)
    if mg_id is None:
      print(row.Symbol)
      mg_id = pymysql_cursor('INSERT INTO Meta.Gene(Symbol) values("%s");' % row.Symbol)
      mg_id = pymysql_cursor('SELECT LAST_INSERT_ID();')
  gene2meta.append([mg_id, row.GeneID, 'Gene'])

Gene2Meta = pd.DataFrame(gene2meta, columns=['MetaID', 'DomainID', 'Class'])



##### Variant
Variant = cursor.execute('SELECT GeneID, VariantID, VariantOrHaplotype From ClinAnn;')
Variant = cursor.fetchall()
Variant = pd.DataFrame(Variant, columns = ["GeneID", "VariantID", "Name"]).drop_duplicates()

rsids = []
variant2meta = []
for index, row in Variant.iterrows():
  var = row.Name
  if var == 'rs58425034':
    var = 'rs12721650'
  elif var == 'rs36056065':
    var = 'rs35854239'
  elif var == 'rs1799735':
    continue
  
  mg_id = Gene2Meta[Gene2Meta.DomainID == row.GeneID].MetaID.to_list()[0]
  if var.startswith('rs'):
    mv_id = pymysql_cursor('SELECT ID FROM Meta.Variant WHERE rsID = "%s";' % var)
  else:
    mv_id = pymysql_cursor('SELECT ID FROM Meta.Variant WHERE GeneID = "%s" AND Name = "%s";' % (mg_id, var))
  # 如果没有匹配到，则补充进去
  if mv_id is None:
    OriginID = pymysql_cursor('SELECT ID FROM Meta.OriginDic WHERE Name = "germline";')
    AssemblyID = pymysql_cursor('SELECT ID FROM Meta.AssemblyDic WHERE Name = "na";')
    ClinicalSigID = pymysql_cursor('SELECT ID FROM Meta.ClinicalSignificanceDic WHERE Name = "drug response";')
    # Symbol
    Symbol = pymysql_cursor('SELECT Symbol FROM Gene WHERE ID = "%s";' % row.GeneID)
    ## 根据是否为rsID分别处理
    # 第一种，一般是单倍型
    if var.startswith('rs') is False and '*' in var:
      TypeID = pymysql_cursor('SELECT ID FROM Meta.VariantTypeDic WHERE Name = "haplotype";')
      mv_id = pymysql_cursor('INSERT INTO Meta.Variant (Symbol, GeneID, Name, TypeID, OriginID, AssemblyID, ClinicalSigID) \
        values ("%s", "%s", "%s", "%s", "%s", "%s", "%s");' % (Symbol, mg_id, var, TypeID, OriginID, AssemblyID, ClinicalSigID))
      mv_id = pymysql_cursor('SELECT LAST_INSERT_ID();')
    # 第二种，rsID类别的需要单独保存下来，然后使用dbSNP的API查询并补充
    elif var.startswith('rs'):
      rsids.append(var)
    # 第三种，是一种及其少见的变异
    else:
      TypeID = pymysql_cursor('SELECT ID FROM Meta.VariantTypeDic WHERE Name = "other";')
      mv_id = pymysql_cursor('INSERT INTO Meta.Variant (Symbol, GeneID, Name, TypeID, OriginID, AssemblyID, ClinicalSigID) \
        values ("%s", "%s", "%s", "%s", "%s", "%s", "%s");' % (Symbol, mg_id, var, TypeID, OriginID, AssemblyID, ClinicalSigID))
      mv_id = pymysql_cursor('SELECT LAST_INSERT_ID();')
  
  # Collect the results
  if type(mv_id) is not list:
    mv_id = [mv_id]
  for mv in mv_id:
    variant2meta.append((mv, row.VariantID, 'Variant'))

print(rsids) # rsids需要人工增加后再继续匹配

Variant2Meta = pd.DataFrame(variant2meta, columns=['MetaID', 'DomainID', 'Class'], dtype='object')



##### Drug
Drug = cursor.execute('SELECT DrugID, Drug From ClinAnn;')
Drug = cursor.fetchall()
Drug = pd.DataFrame(Drug, columns = ["DrugID", "Name"]).drop_duplicates()

drug2meta = []
for index, row in Drug.iterrows():
  ### Use Name
  Name = row.Name; NameTrim = drug_trim(row.Name)
  mdrug_id = map_meta(Name, NameTrim, 'Meta.Drug', 'Meta.DrugSynonym', False)
  ### Cross-Ref
  if mdrug_id is None:
    crossref = pymysql_cursor('SELECT CrossRef FROM Drug WHERE ID = "%s";' % row.DrugID)
    if crossref:
      refs = crossref.strip().replace('"', '').split(',')
      for ref in refs:
        Source = ref.split(':', 1)[0]
        Accession = ref.split(':', 1)[1]
        mdrug_id = pymysql_cursor('SELECT DrugID FROM Meta.DrugExternalLink WHERE Source = "%s" AND Accession = "%s";' % (Source, Accession))
        if mdrug_id:
          break
  ### Use GenericName, TradeName, Synonyms
  if mdrug_id is None:
    Synonyms = []
    generic_name = pymysql_cursor('SELECT GenericName FROM Drug WHERE ID = "%s";' % row.DrugID)
    trade_name = pymysql_cursor('SELECT TradeName FROM Drug WHERE ID = "%s";' % row.DrugID)
    if generic_name:
      Synonyms.append(generic_name)
    if trade_name:
      Synonyms.append(trade_name)
    syn = pymysql_cursor('SELECT Synonym FROM DrugSynonyms WHERE DrugID = "%s";' % row.DrugID)
    if syn:
      if type(syn) != list:
        Synonyms.append(syn)
      else:
        Synonyms.extend(syn)
    # One by one searching
    for Synonym in Synonyms:
      SynonymTrim = drug_trim(Synonym)
      mdrug_id = map_meta(Synonym, SynonymTrim, 'Meta.Drug', 'Meta.DrugSynonym', False)
      if mdrug_id:
        break
  ### Insert
  if mdrug_id is None:
    print("None-map: ", Name, mdrug_id)
    insert_row = pymysql_cursor('INSERT INTO Meta.Drug(Name, NameTrim) values("%s", "%s");' % (Name, NameTrim))
    mdrug_id = pymysql_cursor('SELECT LAST_INSERT_ID();')
    for Synonym in Synonyms:
      SynonymTrim = drug_trim(Synonym)
      insert_row = pymysql_cursor('INSERT INTO Meta.DrugSynonym(Synonym, SynonymTrim, DrugID) values("%s", "%s", "%s");' % (Synonym, SynonymTrim, mdrug_id))
  
  # Collect the results
  if type(mdrug_id) is not list:
    mdrug_id = [mdrug_id]
  for mdrug in mdrug_id:
    drug2meta.append((mdrug, row.DrugID, 'Drug'))


Drug2Meta = pd.DataFrame(drug2meta, columns=['MetaID', 'DomainID', 'Class'])




##### Disease
# 将Type为Not cancer的项过滤掉，仅保留Cancer以及-
Disease = cursor.execute('SELECT ID, Name From Phenotype WHERE Type = "Cancer" OR Type = "-";')
Disease = cursor.fetchall()
Disease = pd.DataFrame(Disease, columns = ["PhenotypeID", "Name"]).drop_duplicates()

disease2meta = []
for index, row in Disease.iterrows():
  ### Use Name
  Name = row.Name; NameTrim = disease_trim(Name)
  mdis_id = map_meta(Name, NameTrim, 'Meta.Disease', 'Meta.DiseaseSynonym', False)
  ### Use Synonym
  if mdis_id is None:
    Synonyms = pymysql_cursor('SELECT Synonym FROM PhenotypeSynonyms WHERE PhenotypeID = "%s";' % row.PhenotypeID)
    if Synonyms:
      if type(Synonyms) != list:
        Synonyms = [Synonyms]
      for Synonym in Synonyms:
        SynonymTrim = disease_trim(Synonym)
        mdis_id = map_meta(Synonym, SynonymTrim, 'Meta.Disease', 'Meta.DiseaseSynonym', False)
        if mdis_id:
          break
  if mdis_id is None:
    print("None-map: ", Name, mdis_id)
    mdis_id = pymysql_cursor('INSERT INTO Meta.Disease(Name, NameTrim) values("%s", "%s");' % (Name, NameTrim))
    mdis_id = pymysql_cursor('SELECT LAST_INSERT_ID();')
  
  # Collect the results
  if type(mdis_id) is not list:
    mdis_id = [mdis_id]
  for mdis in mdis_id:
    disease2meta.append((mdis, row.PhenotypeID, 'Disease'))


Disease2Meta = pd.DataFrame(disease2meta, columns=['MetaID', 'DomainID', 'Class'])

# # Should be Null
# pymysql_cursor('SELECT ID FROM PharmGKB.Phenotype WHERE Type = "Cancer" and MetaDiseaseID is NULL;')

# # Not Cancer
# not_cancer_id = pymysql_cursor('SELECT ID FROM PharmGKB.Phenotype WHERE Type = "Not cancer" and MetaDiseaseID is NULL;')
# md_id = pymysql_cursor('SELECT ID FROM Meta.Disease WHERE Name = "Not cancer";')
# pymysql_cursor('UPDATE Phenotype SET MetaDiseaseID = "%s" WHERE ID IN (%s);' % (md_id, ','.join([str(x) for x in not_cancer_id])))

# # Some clinical annotataions don't have phenotype
# md_id = pymysql_cursor('SELECT ID FROM Meta.Disease WHERE Name is NULL;')
# pymysql_cursor('UPDATE Phenotype SET MetaDiseaseID = "%s" WHERE Phenotype = "";' % md_id)


### Merge the above data
Domain2Meta = pd.concat([Gene2Meta, Variant2Meta, Drug2Meta, Disease2Meta], axis=0)
Domain2Meta.to_csv(TSV['Domain2Meta'], sep='\t', index=0)
# Load data
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Domain2Meta FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (MetaID, DomainID, Class);" % TSV['Domain2Meta'])
logger.info('Domain2Meta records: %i' % row_count)
