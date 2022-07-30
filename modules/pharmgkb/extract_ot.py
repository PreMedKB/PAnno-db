#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract Semantic Network and Therapy Details.
"""

# Import modules
import toml, re
import pandas as pd
from utils.pgi_functions import connect_database, pymysql_cursor, load_logging_cfg
from utils.ot_functions import extract_ot_therapy, get_norm_level, insert_otsnu

cfile = './conf/pharmgkb.toml'
cf = toml.load(cfile, _dict=dict)
global DEFAULT
DEFAULT = cf['DEFAULT']
db, cursor = connect_database(DEFAULT)
logger = load_logging_cfg(DEFAULT)


#########################
## Exact Ontology
#########################
# Fetch PharmGKB treatment data: only consider cancer or '-' phenotype
database = "PharmGKB"
SourceID = pymysql_cursor('SELECT ID FROM OT.SourceDic WHERE Name = "%s";' % database)
ClinAnn = cursor.execute('SELECT * FROM ClinAnn WHERE ID IN (SELECT ClinAnnID FROM ClinAnnHasPhenotype WHERE PhenotypeID IN (SELECT ID FROM Phenotype WHERE Type = "Cancer" OR Type = "-"));')
ClinAnn = cursor.fetchall()
colnames = pymysql_cursor('SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE table_name="ClinAnn" AND table_schema="PharmGKB";')
ClinAnn = pd.DataFrame(ClinAnn, columns=colnames)


# Start to prase the OT relationships!
for index, row in ClinAnn.iterrows():
  
  print("%s/%s" % (index+1, ClinAnn.shape[0]))
  EntryID = row.ID
  OriginID = pymysql_cursor('SELECT ID FROM OT.OriginDic WHERE Name = "Germline";')
  
  # EvidenceLevelIDs
  raw_level = row.EvidenceLevel
  RawLevel = "%s_%s" % (database, raw_level)
  NormLevel = get_norm_level(RawLevel)
  RawLevelID = pymysql_cursor('SELECT ID FROM OT.EvidenceLevelDic WHERE Name = "%s";' % RawLevel)
  NormLevelID = pymysql_cursor('SELECT ID FROM OT.EvidenceLevelDic WHERE Name = "%s";' % NormLevel)
  ## 更新NormLevelID
  # pymysql_cursor('UPDATE OT.TherapyDetail SET NormLevelID = "%s" WHERE SourceID = "%s" AND EntryID = "%s";' % (NormLevelID, SourceID, EntryID))

  # Support
  if raw_level != "4":
    Support = 1
  else:
    Support = 0
  
  Alteration = row.VariantOrHaplotype
  Allele = row.GenotypeOrAllele
  TumorType = row.Phenotypes
  Therapy = row.Drug
  
  # Response: 将反应类别和具体变化拼接起来
  PhenotypeCategory = pymysql_cursor('SELECT Name FROM PhenotypeCategoryDic WHERE ID = "%s";' % row.PhenotypeCategoryID)
  if re.search('Low|Decreased', row.Function):
    Function = "decreased"
  elif re.search('Moderate', row.Function):
    Function = "moderate"
  elif re.search('High|Increased', row.Function):
    Function = "increased"
  elif re.search('Uncertain|Unknown', row.Function):
    Function = "uncertain"
  elif re.search('Normal|No', row.Function):
    Function = "insusceptible"
  
  Response = "%s, %s" % (PhenotypeCategory, Function)
  URL = row.URL
  
  #### OT.TherapyDetail ####
  DetailID = pymysql_cursor('SELECT ID FROM OT.TherapyDetail WHERE SourceID = "%s" AND EntryID = "%s" AND Support = "%s" AND Alteration = "%s" AND Allele = "%s" AND TumorType = "%s" AND Therapy = "%s" AND Response = "%s" AND URL = "%s" AND OriginID = "%s" AND RawLevelID = "%s" AND NormLevelID = "%s";' % (SourceID, EntryID, Support, Alteration, Allele, TumorType, Therapy, Response, URL, OriginID, RawLevelID, NormLevelID))
  if DetailID is None:
    ins_row = pymysql_cursor('INSERT INTO OT.TherapyDetail (SourceID, EntryID, Support, Alteration, Allele, TumorType, Therapy, Response, URL, OriginID, RawLevelID, NormLevelID) VALUES \
      ("%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s");' % \
        (SourceID, EntryID, Support, Alteration, Allele, TumorType, Therapy, Response, URL, OriginID, RawLevelID, NormLevelID))
    DetailID = pymysql_cursor('SELECT LAST_INSERT_ID();')
  
  #### OT.ClinicalAnnotation ####
  Annotation = row.Annotation.strip('"|\'').replace('"', "'")
  if Annotation != '':
    Annotation = "Phenotype Function: %s" % Annotation
    ins_row = pymysql_cursor('INSERT INTO OT.ClinicalAnnotation (DetailID, Annotation) SELECT "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.ClinicalAnnotation WHERE DetailID = "%s" and Annotation = "%s");' % (DetailID, Annotation, DetailID, Annotation))
  
  #### OT.Guideline #### 
  DosingGuideline = cursor.execute('SELECT Source, Annotation FROM ClinDosingGuideline WHERE RelatedGeneID = "%s" AND RelatedDrugID = "%s";' % (row.GeneID, row.DrugID))
  DosingGuideline = cursor.fetchall()
  if DosingGuideline:
    for ele in DosingGuideline:
      Guideline = "%s: %s" % (ele[0], ele[1].strip('"').replace('""', "'").replace('"', "'"))
      ins_row = pymysql_cursor('INSERT INTO OT.Guideline (DetailID, Guideline) SELECT "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.Guideline WHERE DetailID = "%s" and Guideline = "%s");' % (DetailID, Guideline, DetailID, Guideline))
  
  #### OT.RefPubMed ####
  pubmeds = pymysql_cursor('SELECT PMID FROM ClinAnnEvidence WHERE CAID = "%s";' % row.CAID)
  if type(pubmeds) is not list:
    pubmeds = [pubmeds]
  for pmid in pubmeds:
    if pmid != '':
      ins_row = pymysql_cursor('INSERT INTO OT.RefPubMed(DetailID, PMID) SELECT "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.RefPubMed WHERE DetailID = "%s" and PMID = "%s");' % (DetailID, pmid, DetailID, pmid))
  
  #### OT.Therapy ####
  # PharmGKB比较简单，不涉及组合的问题
  VariantSetID = None
  DrugSetID = None
  mvar_id = pymysql_cursor('SELECT DISTINCT MetaID FROM Domain2Meta WHERE Class = "Variant" AND DomainID = "%s";' % row.VariantID)
  mgene_id = pymysql_cursor('SELECT DISTINCT MetaID FROM Domain2Meta WHERE Class = "Gene" AND DomainID = "%s";' % row.GeneID)
  mdrug_id = pymysql_cursor('SELECT DISTINCT MetaID FROM Domain2Meta WHERE Class = "Drug" AND DomainID = "%s";' % row.DrugID)
  mdis_id = pymysql_cursor('SELECT DISTINCT MetaID FROM Domain2Meta WHERE Class = "Disease" AND DomainID IN (SELECT PhenotypeID FROM ClinAnnHasPhenotype WHERE ClinAnnID = "%s");' % row.ID)
  
  # 不应该有None, 否则说明前面的步骤出错了
  # 2021.12.14有三条没有匹配上，下次再看：30193 9423 None; [211, 3886, 4243, 4363, 9064, 17893] 16945 None; [211, 3886, 4243, 4363, 9064, 17893] 11769 None
  if mdis_id is None or mdrug_id is None or mvar_id is None:
    print(mdis_id, mdrug_id, mvar_id)
    continue
  
  if type(mvar_id) != list:
    mvar_id = [mvar_id]
  if type(mdrug_id) != list:
    mdrug_id = [mdrug_id]
  if type(mdis_id) != list:
    mdis_id = [mdis_id]
  
  # Inserting
  for v_id in mvar_id:
    for dis_id in mdis_id:
      for dr_id in mdrug_id:
        if VariantSetID:
          VariantID = None
        else:
          VariantID = v_id
        if DrugSetID:
          DrugID = None
        else:
          DrugID = dr_id
        DiseaseID = dis_id
        extract_ot_therapy(VariantID, VariantSetID, DrugID, DrugSetID, DiseaseID, DetailID)
  
  #### OT.SemanticNetwork ####
  direction_id = pymysql_cursor("SELECT ID FROM OT.DirectionDic WHERE Name = 'left2right';")
  # Variant-Disease
  if VariantSetID:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'VariantSet';")
  else:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Variant';")
  class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Disease';")
  relations = ["cause"]
  for v_id in mvar_id:
    for dis_id in mdis_id:
      id1 = v_id
      id2 = dis_id
      insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID)
  
  # Variant-Gene
  if VariantSetID:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'VariantSet';")
  else:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Variant';")
  class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Gene';")
  relations = ["locate in and effect"]
  for v_id in mvar_id:
    id1 = v_id
    id2 = mgene_id
    insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID)
  
  # Gene-Disease
  class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Gene';")
  class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Disease';")
  relations = ["oncogene"]
  for dis_id in mdis_id:
    id1 = mgene_id
    id2 = dis_id
    insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID)
  
  # Variant-Drug
  relations = ["germline"]
  if VariantSetID:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'VariantSet';")
  else:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Variant';")
  if DrugSetID:
    class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'DrugSet';")
  else:
    class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Drug';")
  for v_id in mvar_id:
    for dr_id in mdrug_id:
      id1 = v_id
      id2 = dr_id
      insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID)
  
  # Drug-Disease
  if DrugSetID:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'DrugSet';")
  else:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Drug';")
  class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Disease';")
  relations = ["cure"]
  for dr_id in mdrug_id:
    for dis_id in mdis_id:
      id1 = dr_id
      id2 = dis_id
      insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID)
  
  # Drug-Gene
  if DrugSetID:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'DrugSet';")
  else:
    class_id1 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Drug';")
  class_id2 = pymysql_cursor("SELECT ID FROM OT.ClassDic WHERE Name = 'Gene';")
  relations = ['associate']
  for dr_id in mdrug_id:
    id1 = dr_id
    id2 = mgene_id
    insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID)

