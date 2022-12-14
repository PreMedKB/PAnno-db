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
from utils.parse_clinann import infer_response

cfile = './conf/pharmgkb.toml'
cf = toml.load(cfile, _dict=dict)
DEFAULT = cf['DEFAULT']
RAW = cf['RAW']
TSV = cf['TSV']

db, cursor = connect_database(DEFAULT)
logger = load_logging_cfg(DEFAULT)


####### Drug labels ########
import requests, re, json, math
from bs4 import BeautifulSoup
headers={'Cookie': '_ga=GA1.2.1379803965.1662683552; _gid=GA1.2.1271233879.1665052486; _gat=1',
    'Host': 'api.pharmgkb.org', 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/105.0.0.0 Safari/537.36'}
    #'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.0 Safari/605.1.15',

fail = []
drug_label = pd.read_csv(RAW['drug_label_tsv'], sep="\t")
for index, row in drug_label.iterrows():
  print('%s/%s' % (index+1, drug_label.shape[0]))
  PAID = row['PharmGKB ID']
  try:
    url = 'https://api.pharmgkb.org/v1/site/page/overview/{0}?view=max'.format(PAID)
    response = requests.get(url, headers=headers)
    contents = json.loads(response.text)
    file_name = '{0}/{1}.json'.format(RAW['drug_label_path'], re.sub('___', '_', re.sub(' |/', '_', row['Name'])))
    with open(file_name, 'w') as json_file:
      json_str = json.dumps(contents['data'], indent=4)
      json_file.write(json_str)
  except:
    fail.append(PAID)

print(fail)


files = os.listdir(RAW['drug_label_path'])
labels = []
for f in files:
  if ".json" in f:
    Source = f.split("_")[2]
    with open("%s/%s" % (RAW['drug_label_path'], f)) as load_f:
      content = json.load(load_f)
      # Check
      if len(content["objects"]) > 1:
        print(f)
      URL = "https://www.pharmgkb.org/labelAnnotation/{0}".format(content["objects"][0]["id"])
      raw_sum = html.unescape(content["objects"][0]["summaryMarkdown"]["html"]).strip()
      Annotation = re.sub('<em>|</em>|<p>|</p>|<strong>|</strong>|<a href=".*">|</a>|<a rel=".*">', '', raw_sum)
      Annotation = re.sub('\n', '', Annotation)
      # Remove '"", etc.
      Annotation = Annotation.strip('"').replace('""', "'").replace('"', "'")
      # Related Drugs and Genes
      drugs = content["objects"][0]["relatedChemicals"]
      genes = content["objects"][0]["relatedGenes"]
      for drug in drugs:
        DrugID = pymysql_cursor("SELECT ID FROM Drug WHERE PAID ='%s';" % drug["id"])
        Drug = drug['name']
        for gene in genes:
          GeneID = pymysql_cursor("SELECT ID FROM Gene WHERE PAID ='%s';" % gene["id"])
          Gene = gene['symbol']
          labels.append([Source, Annotation, Gene, Drug, GeneID, DrugID, URL])

DrugLabels = pd.DataFrame(labels, columns=["Source", "Summary", "Gene", "Drug", "GeneID", "DrugID", "URL"])


####### Dosing guidelines ########
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
        DrugID = pymysql_cursor("SELECT ID FROM Drug WHERE PAID ='%s';" % drug["id"])
        Drug = drug['name']
        for gene in genes:
          GeneID = pymysql_cursor("SELECT ID FROM Gene WHERE PAID ='%s';" % gene["id"])
          Gene = gene['symbol']
          guidelines.append([Source, Annotation, Gene, Drug, GeneID, DrugID, URL])

ClinDosingGuideline = pd.DataFrame(guidelines, columns=["Source", "Summary", "Gene", "Drug", "GeneID", "DrugID", "URL"])


Guideline = pd.concat([DrugLabels, ClinDosingGuideline], axis = 0)
table_p = TSV["Guideline"]
Guideline.to_csv(table_p, sep="\t", header=0, index=0)
row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE Guideline FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Source, Annotation, Gene, Drug, GeneID, DrugID, URL);" % table_p)
row_count