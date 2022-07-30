#!/usr/bin/env python3
# -*- coding: utf-8 -*-

' Trim of Drug & Disease Names Used in PreMedKB. '

__author__ = 'Yaqing Liu'


import inflect
import re
import nltk
from nltk import SnowballStemmer

###################
## Drug
###################
# Load reference data
drug_ignored_all = [x.strip() for x in open("./utils/assets/drugname_ignored_word_all.txt").readlines()]
drug_ignored = [x.strip() for x in open("./utils/assets/drugname_ignored_word.txt").readlines()]

def drug_trim(drug_name): 
  p = inflect.engine()
  
  # Change to lower case
  name = drug_name.lower()
  
  # Convert plurals to singulars
  singular_name = p.singular_noun(name)
  if singular_name is False:
    singular_name = name
  
  # Treatment of nouns with brackets ()
  singular_name = singular_name.replace(' (', '(')
  within_brackets = re.findall(r'[(](.*?)[)]', singular_name)
  if within_brackets:
    # Some name is like '(DMSO)', in this case, we should not remove the brackets and content with the brackets.
    complete_bra = "(" + within_brackets[0] + ")"
    if complete_bra == singular_name:
      name_rmbracket = within_brackets[0]
    else:
      name_rmbracket = singular_name
      for bra in within_brackets:
        complete_bra = "(" + bra + ")"
        if ";" in complete_bra:
          complete_bra_rp = re.sub("\(|\)", ";", complete_bra) 
          name_rmbracket = name_rmbracket.replace(complete_bra, complete_bra_rp)
        else:
          name_rmbracket = name_rmbracket.replace(complete_bra, "")
  else:
    name_rmbracket = singular_name
  
  # Remove extra spaces. ' +' means have continuous spaces
  name_rmbracket = re.sub(r" +", " ", name_rmbracket).strip()
  
  # drugname_ignored_word_all
  # The list of words that can be removed after matching in any position
  name_list = re.split(r" |;", name_rmbracket)
  if len(name_list) >= 2:
    name_stop = []
    for n in name_list:
      name_stop.append(re.sub("|".join(drug_ignored_all), "", n))
  else:
    name_stop = name_list
  
  # drugname_ignored_word: >= 2 words
  # Match and remove as a whole in word form
  name_ignored = [name_stop[0]]
  if len(name_stop) >= 2:
    for i in range(1, len(name_stop)):
      if name_stop[i] not in drug_ignored:
        name_ignored.append(name_stop[i])
  name_ignored = " ".join(name_ignored)
  
  # Only keep number and alphabet, replace punctuations with space
  name_convert = re.sub(r"[^a-z0-9 ]", " ", name_ignored)

  # Separate numbers and letters by spaces
  wd = re.findall(r'([a-z]+)([0-9]+)', name_convert)
  if wd:
    for i in wd:
      a = ''.join(i)
      b = ' '.join(i)
      name_convert = name_convert.replace(a,b)
  dw = re.findall(r'([0-9]+)([a-z]+)', name_convert)
  if dw:
    for i in dw:
      a = ''.join(i)
      b = ' '.join(i)
      name_convert = name_convert.replace(a,b)
  
  # Remove extra spaces. ' +' means have continuous spaces
  name_trim = re.sub(r" +", " ",name_convert).strip()
  
  # Sort by word
  name_trim = " ".join(sorted(name_trim.split(" ")))
  return(name_trim)




###################
## Disease
###################
# Load reference data
disease_stoplist = [x.strip().lower() for x in open("./utils/assets/disease_stoplist_stem.txt").readlines()]

replacerule_general = {}
for i in [x.strip() for x in open("./utils/assets/disease_replacerule_general.txt").readlines()]:
  i = i.strip().split('|')
  replacerule_general[i[0]] = i[1]

replacerule_cancer = {}
for i in [x.strip() for x in open("./utils/assets/disease_replacerule_cancer.txt").readlines()]:
  i = i.strip().split('|')
  replacerule_cancer[i[0]] = i[1]

def disease_trim(disease_name):
  SnowballStemmer.languages
  stemmer = SnowballStemmer("english")
  stemmer.stem(u"Autobahnen")
  ps = nltk.stem.snowball.PortugueseStemmer()
  
  # Change to lower case
  name = disease_name.lower()
  
  # Convert characters other than [^a-zA-Z0-9_] to spaces
  name_convert = re.sub(r"[^a-z0-9 ]"," ",name)
  
  # Separate numbers and letters by spaces
  wd = re.findall(r'([a-z]+)([0-9]+)', name_convert)
  if wd:
    for i in wd:
      a = ''.join(i)
      b = ' '.join(i)
      name_convert = name_convert.replace(a,b)
  dw = re.findall(r'([0-9]+)([a-z]+)', name_convert)
  if dw:
    for i in dw:
      a = ''.join(i)
      b = ' '.join(i)
      name_convert = name_convert.replace(a,b)
  
  # Replace the term used to represent 'cancer'
  # Note that: we should split the name into list, otherwise the replacement will occur like 'adenocarcinoma' -> 'adenocancer', while 'cancer' is expected.
  name_list = name_convert.split(" ")
  name_list_cancer = [replacerule_cancer[rc] if rc in replacerule_cancer else rc for rc in name_list]
  
  # Replace the general items, such as i, ii, ..., alpha, beta
  # Note that: we should split the name into list, otherwise the replacement will occur within a word which will effect the following 'ps.stem' step.
  name_list_general = [replacerule_general[rg] if rg in replacerule_general else rg for rg in name_list_cancer]

  # Sort by word
  name_sort = [j for j in name_list_general if j != '']
  name_sort.sort()
  
  # Word-by-word stem
  name_stem = []
  for word in name_sort:
    trim_word = ps.stem(word)
    name_stem.append(trim_word)
  
  # Determine if the stemmed word is in the stop list
  name_trim_list = []
  for ns in name_stem:
    if ns not in disease_stoplist:
      name_trim_list.append(ns)
  name_trim = ' '.join(name_trim_list)
  
  return(name_trim)