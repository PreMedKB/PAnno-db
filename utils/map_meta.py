### Function: Mapping to MetaDisease/Drug by names or synonyms
### 完整版可见utils.map_meta
import re

def map_meta(Term, TermTrim, Table_Name, Table_Synonym, Insert=False):
  
  # Clean Term
  Term = Term.strip('"').replace('""', "'").replace('"', "'")
  # Check
  ID_LIST = []
  ELEMENT_LIST = ['Gene', 'Variant', 'Disease', 'Drug']
  for element in ELEMENT_LIST:
    element1 = re.findall(element.lower(), Table_Name.lower())
    element2 = re.findall(element.lower(), Table_Synonym.lower())
    if element1 != []:
      if element1 != element2:
        print('Error: %s does not match with the %s! Please check.' % (Table_Name, Table_Synonym))
      else:
        ID_Name = element.capitalize() + 'ID'
        ID_LIST.append(ID_Name)
  if len(set(ID_LIST)) != 1:
    print('Error: %s conflicts with the elements in Meta! Please check.' % (Table_Name))
  else:
    SynonymID = ID_LIST[0]
  
  # Exact Match
  md_id = pymysql_cursor('SELECT ID FROM %s WHERE Name = "%s";' % (Table_Name, Term))
  if md_id is None:
    md_id = pymysql_cursor('SELECT %s FROM %s WHERE Synonym = "%s";' % (SynonymID, Table_Synonym, Term))
  if md_id is None:
    md_id = pymysql_cursor('SELECT ID FROM %s WHERE NameTrim = "%s";' % (Table_Name, TermTrim))
  if md_id is None:
    md_id = pymysql_cursor('SELECT %s FROM %s WHERE SynonymTrim = "%s";' % (SynonymID, Table_Synonym, TermTrim))
  
  # Fuzzy Match
  if SynonymID == 'DrugID':
    min_length = 5
  elif SynonymID == 'DiseaseID':
    min_length = 3
  if len(TermTrim.split()) > min_length:
    term = '%s%s%s' % ('%', Term, '%')
    termtrim = '%s%s%s' % ('%', TermTrim, '%')
    if md_id is None:
      md_id = pymysql_cursor('SELECT ID FROM %s WHERE Name LIKE "%s";' % (Table_Name, term))
    if md_id is None:
      md_id = pymysql_cursor('SELECT %s FROM %s WHERE Synonym LIKE "%s";' % (SynonymID, Table_Synonym, term))
    if md_id is None:
      md_id = pymysql_cursor('SELECT ID FROM %s WHERE NameTrim LIKE "%s";' % (Table_Name, termtrim))
    if md_id is None:
      md_id = pymysql_cursor('SELECT %s FROM %s WHERE SynonymTrim LIKE "%s";' % (SynonymID, Table_Synonym, termtrim))
  
  return(md_id)
