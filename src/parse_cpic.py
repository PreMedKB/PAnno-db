import os, toml
import pandas as pd
from utils.pgi_functions import connect_database, pymysql_cursor, load_logging_cfg

cfile = './conf/panno.toml'
cf = toml.load(cfile, _dict=dict)
DEFAULT = cf['DEFAULT']
RAW = cf['RAW']
TSV = cf['TSV']

db, cursor = connect_database(DEFAULT)
logger = load_logging_cfg(DEFAULT)

gene_list = ["ABCG2", "CACNA1S", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "CFTR", "DPYD", "G6PD", "HLA-A", "HLA-B", "IFNL3", "MT-RNR1", "NUDT15", "RYR1", "SLCO1B1", "TPMT", "UGT1A1", "VKORC1", "CYP3A4", "CYP2C8"]

# definition = './data/cpic/allele_definition'
# functionality = './data/cpic/allele_functionality'
# diplotype_phenotype = './data/cpic/diplotype_phenotype'
# frequency = './data/cpic/frequency'
# example_cds = './data/cpic/example_cds'

functionality = pd.DataFrame()
diplotype_phenotype = pd.DataFrame()
example_cds = pd.DataFrame()
for gene in gene_list:
  # f1p = './data/cpic/allele_definition/{0}_allele_definition_table.xlsx'.format(gene)
  # if os.path.isfile(f1p):
  #   definition = pd.read_excel(f1p, skiprows=1)
  
  f2p = './data/cpic/allele_functionality/{0}_allele_functionality_reference.xlsx'.format(gene)
  if os.path.isfile(f2p):
    df2 = pd.read_excel(f2p, skiprows=1, usecols = [0,1,3,5])
    df2.insert(0, 'Gene', gene)
    df2.columns = ['Gene', 'Allele', 'ActivityScore', 'Function', 'PMID']
    functionality = pd.concat([functionality, df2], axis = 0).reset_index(drop = True).fillna('')
  
  f3p = './data/cpic/diplotype_phenotype/{0}_Diplotype_Phenotype_Table.xlsx'.format(gene)
  if os.path.isfile(f3p):
    df3 = pd.read_excel(f3p, usecols = [0,1,2,3])
    df3.insert(0, 'Gene', gene)
    # DPYD Telti/Kobe
    df3 = pd.concat([df3, df3.iloc[:,1].str.split('/', expand=True, n=1)], axis = 1)
    if gene == 'MT-RNR1':
      df3.insert(6, 'Allele2', '')
    df3.columns = ['Gene', 'Diplotype', 'ActivityScore', 'Phenotype', 'EHR', 'Allele1', 'Allele2']
    diplotype_phenotype = pd.concat([diplotype_phenotype, df3], axis = 0).reset_index(drop = True).fillna('')
  
  f4p = './data/cpic/example_cds/{0}_CDS.xlsx'.format(gene)
  if os.path.isfile(f4p):
    df4 = pd.read_excel(f4p, skiprows=1, usecols = [0,1,2,3])
    df4.insert(0, 'Gene', gene)
    df4.columns = ['Gene', 'Phenotype', 'ActivityScore', 'EHR', 'Consultation']
    example_cds = pd.concat([example_cds, df4], axis = 0).reset_index(drop = True).fillna('')


##### Table EHRConsultation
table_p = TSV["EHRConsultation"]
example_cds.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE EHRConsultation FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Gene, Phenotype, ActivityScore, EHR, Consultation);" % table_p)
row_count


##### Table DiplotypePhenotype
diplotype_phenotype.insert(diplotype_phenotype.shape[1], 'ConsultationID', '')
for index, row in diplotype_phenotype.iterrows():
  ConsultationID = pymysql_cursor("SELECT ID FROM EHRConsultation WHERE Gene = '%s' AND Phenotype = '%s' AND ActivityScore = '%s' AND EHR = '%s';" % (row.Gene, row.Phenotype, row.ActivityScore, row.EHR))
  if ConsultationID is None:
    print(index)
  elif type(ConsultationID) is list:
    print(index)
  row.ConsultationID = ConsultationID
  diplotype_phenotype.iloc[index] = row

table_p = TSV["DiplotypePhenotype"]
diplotype_phenotype = diplotype_phenotype[['Gene', 'Allele1', 'Allele2', 'ActivityScore', 'Phenotype', 'EHR', 'ConsultationID']]
diplotype_phenotype.to_csv(table_p, sep="\t", index=0, quoting=3)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE DiplotypePhenotype FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Gene, Allele1, Allele2, ActivityScore, Phenotype, EHR, ConsultationID);" % table_p)
row_count


##### Table AlleleFunctionality
# Save and manually modified the data
functionality.insert(functionality.shape[1], 'Variant', '')
functionality.insert(functionality.shape[1], 'AlleleManual', '')
functionality.insert(functionality.shape[1], 'FunctionManual', '')
# functionality.to_csv("./data/cpic/allele_functionality.txt", sep="\t", index=0, encoding='utf-8')

functionality = pd.read_csv("./data/cpic/allele_functionality.txt", sep="\t")
table_p = TSV["AlleleFunctionality"]
functionality.to_csv(table_p, sep="\t", index=0, encoding='utf-8', quoting=3)
# Manual processing ......
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE AlleleFunctionality FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' IGNORE 1 LINES (Gene, Allele, ActivityScore, Function, PMID, Variant, AlleleManual, FunctionManual);" % table_p)
row_count
