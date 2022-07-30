# Import modules
from utils.pgi_functions import pymysql_cursor


# Function: insert new therapy into tables
def extract_ot_therapy(VariantID, VariantSetID, DrugID, DrugSetID, DiseaseID, DetailID):
  # Single Variant
  if VariantSetID is None:
    if DrugSetID is None:
      ins_row = pymysql_cursor('INSERT INTO OT.Therapy(VariantID, DrugID, DiseaseID) SELECT "%s", "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.Therapy WHERE VariantID="%s" and DrugID="%s" and DiseaseID="%s");' % (VariantID, DrugID, DiseaseID, VariantID, DrugID, DiseaseID))
      TherapyID = pymysql_cursor('SELECT * FROM OT.Therapy WHERE VariantID="%s" and DrugID="%s" and DiseaseID="%s";' % (VariantID, DrugID, DiseaseID))
    elif DrugSetID:
      ins_row = pymysql_cursor('INSERT INTO OT.Therapy(VariantID, DrugSetID, DiseaseID) SELECT "%s", "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.Therapy WHERE \
          VariantID="%s" and DrugSetID="%s" and DiseaseID="%s");' % (VariantID, DrugSetID, DiseaseID, VariantID, DrugSetID, DiseaseID))
      TherapyID = pymysql_cursor('SELECT * FROM OT.Therapy WHERE VariantID="%s" and DrugSetID="%s" and DiseaseID="%s";' % (VariantID, DrugSetID, DiseaseID))
    
    # MetaGeneID(s)
    mg_ids = pymysql_cursor('SELECT GeneID FROM MetaLite.Variant WHERE ID = "%s";' % VariantID)
    
  # Variant Set
  elif VariantSetID:
    if DrugSetID is None:
      ins_row = pymysql_cursor('INSERT INTO OT.Therapy(VariantSetID, DrugID, DiseaseID) SELECT "%s", "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.Therapy WHERE VariantSetID="%s" and DrugID="%s" and DiseaseID="%s");' % (VariantSetID, DrugID, DiseaseID, VariantSetID, DrugID, DiseaseID))
      TherapyID = pymysql_cursor('SELECT * FROM OT.Therapy WHERE VariantSetID="%s" and DrugID="%s" and DiseaseID="%s";' % (VariantSetID, DrugID, DiseaseID))
    elif DrugSetID:
      ins_row = pymysql_cursor('INSERT INTO OT.Therapy(VariantSetID, DrugSetID, DiseaseID) SELECT "%s", "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.Therapy WHERE \
          VariantSetID="%s" and DrugSetID="%s" and DiseaseID="%s");' % (VariantSetID, DrugSetID, DiseaseID, VariantSetID, DrugSetID, DiseaseID))
      TherapyID = pymysql_cursor('SELECT * FROM OT.Therapy WHERE VariantSetID="%s" and DrugSetID="%s" and DiseaseID="%s";' % (VariantSetID, DrugSetID, DiseaseID))
    
    # MetaGeneIDs
    mg_ids = pymysql_cursor('SELECT GeneID FROM MetaLite.Variant WHERE ID IN (SELECT MetaID FROM MetaLite.SubsetHasElement WHERE SubsetID IN (SELECT SubSetID FROM MetaLite.SetHasSubset WHERE SetID = "%s"));' % VariantSetID)
  
  #### OT.TherapyTarget ####
  if mg_ids is None:
    print(VariantID)
  if type(mg_ids) != list:
    mg_ids = [mg_ids]
  for mg_id in mg_ids:
    ins_row = pymysql_cursor('INSERT INTO OT.TherapyTarget (GeneID, TherapyID) VALUES ("%s", "%s");' % (mg_id, TherapyID))
  
  #### OT.TherapyHasDetail ####  
  ins_row = pymysql_cursor('INSERT INTO OT.TherapyHasDetail (TherapyID, DetailID) SELECT "%s", "%s" FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.TherapyHasDetail WHERE TherapyID="%s" and DetailID="%s");' % (TherapyID, DetailID, TherapyID, DetailID))



# Function: insert new ot relationship into tables
def insert_otsnu(direction_id, class_id1, class_id2, id1, id2, relations, DetailID):
  # OT.SemanticNetwork
  ins_row = pymysql_cursor("INSERT INTO OT.SemanticNetwork(ID1, ClassID1, ID2, ClassID2, DirectionID) SELECT '%s', '%s', '%s', '%s', '%s' FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.SemanticNetwork WHERE ID1 = '%s' and ClassID1 = '%s' and ID2 = '%s' and ClassID2 = '%s' and DirectionID = '%s');" % (id1, class_id1, id2, class_id2, direction_id, id1, class_id1, id2, class_id2, direction_id))
  SNID = pymysql_cursor("SELECT ID FROM OT.SemanticNetwork WHERE ID1 = '%s' and ClassID1 = '%s' and ID2 = '%s' and ClassID2 = '%s' and DirectionID = '%s';" % (id1, class_id1, id2, class_id2, direction_id))
  # OT.SNHasRelation
  for r in relations:
    relation_id = pymysql_cursor("SELECT ID FROM OT.RelationDic WHERE Name = '%s';" % r)
    ins_row = pymysql_cursor("INSERT INTO OT.SNHasRelation(SNID, RelationID) SELECT '%s', '%s' FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.SNHasRelation WHERE SNID = '%s' and RelationID = '%s');" % (SNID, relation_id, SNID, relation_id))
  # OT.SNHasDetail
  ins_row = pymysql_cursor("INSERT INTO OT.SNHasDetail(SNID, DetailID) SELECT '%s', '%s' FROM DUAL WHERE NOT EXISTS (SELECT * FROM OT.SNHasDetail WHERE SNID = '%s' and DetailID = '%s');" % (SNID, DetailID, SNID, DetailID))



# Fuction: normalize the evidence level
# mapping to AMP/ASCO/CAP Consensus Recommendation
def get_norm_level(raw_level):
  source = raw_level.split('_')[0]
  level = raw_level.split('_')[1]
  # OncoKB
  if source == "OncoKB":
    map_rule = {"1": "A", "2": "A", "R1": "A", "3A": "B", "3B": "C", "4": "D", "R2": "D"}
  # CIViC
  elif source == "CIViC":
    map_rule = {"A": "A", "B": "B", "C": "C", "D": "D", "E": "D"}
  # My Cancer Genome
  elif source == "MCG":
    map_rule = {"Clinical Guidelines": "A", "MCG": "B"}
  # CGI
  elif source == "CGI":
    map_rule = {"Clinical Guidelines": "A", "Late Trials": "B", "Early Trials": "C", "Case Reports": "C", "Preclinical Data": "D", "-": "-"}
    # map_rule = {"Clinical guidelines": "A", "FDA guidelines": "A", "NCCN guidelines": "A", "European LeukemiaNet guidelines": "A", "Late trials": "B", "Early trials": "C", "Case report": "C", "Pre-clinical": "D", "Clinical trials": "-", "Late trials,Pre-clinical": "-"}
  # PharmGKB
  elif source == "PharmGKB":
    map_rule = {"1A": "A", "1B": "A", "2A": "B", "2B": "B", "3": "C",  "4": "-"}
  
  norm_level = map_rule[level]
  if norm_level != "-":
    output = "PreMedKB_%s" % norm_level
  else:
    output = "-"
  
  return(output)

