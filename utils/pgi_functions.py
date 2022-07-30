#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
The functions used in PGIKB.
"""

# import modules
import sys
import toml
import pymysql
import re
import logging.config
from warnings import filterwarnings


# Config
def db_config(cfile):
    global DEFAULT, JSON, TSV, SCHEMA
    cf = toml.load(cfile, _dict=dict)
    DEFAULT = cf['DEFAULT']
    JSON = cf['JSON']
    TSV = cf['TSV']
    SCHEMA = cf['SCHEMA']

# Connect database
def connect_database(DEFAULT):
    # Connect mysql
    global db, cursor
    try:
        db = pymysql.connect(host=DEFAULT['host'], user=DEFAULT['user'], port=DEFAULT['port'],
                             passwd=DEFAULT['passwd'], charset=DEFAULT['charset'], local_infile=DEFAULT['local_infile'])
        # Create a cursor object that can execute SQL statements
        cursor = db.cursor()
        cursor.execute("USE `%s`;" % DEFAULT['database_name'])
    except pymysql.Error as e:
        db.rollback()
        logger.error(e)
        # Define "3" as the db connecting error
        sys.exit(3)
    except pymysql.Warning as w:
        db.rollback()
        logger.warning(w)

    return(db, cursor)


# Execute SQL command; Capture errors and warnings
def pymysql_cursor(sql):
    # try... except... can only capture errors
    filterwarnings("error", category=pymysql.Warning)
    try:
        # cursor.execute(sql) will automatically return the number of affected rows
        effected_rows = cursor.execute(sql)
        db.commit()
    except pymysql.Error as e:
        db.rollback()
        logger.error(e)
    except pymysql.Warning as w:
        db.rollback()
        logger.warning(w)
    else:
        if re.match('^SELECT', sql):
            # cursor.fetchone() return tuple; cursor.execute(sql) returns the number of rows affected
            # one query result: ((1,)); multiple results: ((1,),(2,))
            result = cursor.fetchall()
            if len(result) == 1:
                single_result = result[0][0]
                return single_result
            elif len(result) > 1:
                multiple_results = []
                for ele in result:
                    multiple_results.append(ele[0])
                return multiple_results
        elif re.match('^LOAD DATA LOCAL INFILE', sql):
            return cursor.rowcount


# Logger
def load_logging_cfg(DEFAULT):
    # Define three log output formats
    standard_format = '%(asctime)-15s  [%(threadName)s:%(thread)d, task_id:%(name)s, %(filename)s:%(lineno)d]' \
        '\n[%(levelname)s]  %(message)s\n'
    simple_format = '[%(levelname)s]  %(asctime)-15s\n%(message)s\n'
    #id_simple_format = '[%(levelname)s][%(asctime)s]%(message)s'
    # Full path of log file
    logfile = DEFAULT['logfile']
    # Configuration of the logger
    LOGGING_DIC = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'standard': {
                'format': standard_format
            },
            'simple': {
                'format': simple_format
            },
        },
        'filters': {},
        'handlers': {
            # Print logs to terminal
            'console': {
                'level': 'DEBUG',
                'class': 'logging.StreamHandler',
                'formatter': 'simple',
            },
            # Print logs into file, collect logs above 'INFO'
            'default': {
                'level': 'DEBUG',
                'class': 'logging.handlers.RotatingFileHandler',
                'formatter': 'standard',
                'filename': logfile,
                'maxBytes': 1024*1024*5,
                'backupCount': 5,
                'encoding': 'utf-8',
            },
        },
        'loggers': {
            # Using logging.getLogger(__name__) to get the configuration
            '': {
                'handlers': ['default', 'console'],
                'level': 'DEBUG',
                'propagate': True,
            },
        },
    }
    logging.config.dictConfig(LOGGING_DIC)
    global logger
    logger = logging.getLogger(__name__)
    return(logger)


# Map meta
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
