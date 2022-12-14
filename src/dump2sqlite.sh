#!/bin/bash

# pip install mysql-to-sqlite3
# mysql2sqlite --help
cd ./data/panno/
rm pgx_kb.sqlite3
mysql2sqlite -d PAnno -u root -f pgx_kb.sqlite3 -t ClinAnn GuidelineMerge GuidelineRule DiplotypePhenotype # AlleleFunctionality EHRConsultation