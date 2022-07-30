#!/bin/bash

# pip install mysql-to-sqlite3
# mysql2sqlite --help
rm pgx_kb.sqlite3
mysql2sqlite -d PharmGKB -u root -f pgx_kb.sqlite3 -t ClinAnn PhenotypeCategoryDic ClinDosingGuideline Drug