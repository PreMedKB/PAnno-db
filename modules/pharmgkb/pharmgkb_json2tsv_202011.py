# coding=utf-8

# Import modules
import sys
import getopt
import toml
import json
import pymysql
import re
import logging.config
from warnings import filterwarnings

# Main
def main(argv):
    cfile = ''
    try:
        opts, args = getopt.getopt(
            argv, "hc:", ["cfile="])
    except getopt.GetoptError:
        print('json2tsv.py -c <congfigure file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(
                'json2tsv.py -c <congfigure file>')
            sys.exit()
        elif opt in ("-c", "--cfile"):
            cfile = arg
    
    config(cfile)
    load_logging_cfg()
    connect_database() 
    json2tsv()

### Config information
def config(cfile):
    # Parse config.toml
    #cfile = '/Users/yaqing/Desktop/PGIKB/PharmGKB/config.toml'
    cfile = '/PGIKB/PharmGKB/config.toml'
    cf = toml.load(cfile, _dict=dict)
    # Paths to json/tsv files and MySQL schema
    global DEFAULT, JSON, TSV, PREPRO, SCHEMA
    DEFAULT = cf['DEFAULT']
    JSON = cf['JSON']
    TSV = cf['TSV']
    PREPRO = cf['PREPRO']
    SCHEMA = cf['SCHEMA']

def connect_database(): 
    # Connect mysql
    global db, cursor
    try:
        db = pymysql.connect(host = DEFAULT['host'], user = DEFAULT['user'], port = DEFAULT['port'], charset = DEFAULT['charset'], local_infile = DEFAULT['local_infile'])
        # Create a cursor object that can execute SQL statements
        cursor = db.cursor()
        cursor.execute("USE %s;" % DEFAULT['database_name'])
    except pymysql.Error as e:
        db.rollback()
        logger.error(e)
        # Define "3" as the db connecting error
        sys.exit(3)
    except pymysql.Warning as w:
        db.rollback()
        logger.warning(w)

### Execute SQL command; Capture errors and warnings
def pymysql_cursor(sql):
    # try... except... can only capture errors
    filterwarnings("error", category = pymysql.Warning)
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
        
### Json to tsv
def json2tsv():
    # Read json files
    clinical_ann_metadata = JSON['clinical_ann_metadata']
    file = open(clinical_ann_metadata).read()
    temp = json.loads(file)
    clinical_ann_metadata = temp['data']

    clinical_ann = JSON['clinical_ann']
    file = open(clinical_ann).read()
    temp = json.loads(file)
    clinical_ann = temp['data']

    study_parameters = JSON['study_parameters']
    file = open(study_parameters).read()
    temp = json.loads(file)
    study_parameters = temp['data']

    var_drug_ann = JSON['var_drug_ann']
    file = open(var_drug_ann).read()
    temp = json.loads(file)
    var_drug_ann = temp['data']

    var_fa_ann = JSON['var_fa_ann']
    file = open(var_fa_ann).read()
    temp = json.loads(file)
    var_fa_ann = temp['data']

    var_pheno_ann = JSON['var_pheno_ann']
    file = open(var_pheno_ann).read()
    temp = json.loads(file)
    var_pheno_ann = temp['data']

    var_ann = var_drug_ann + var_fa_ann + var_pheno_ann

    # Create abnormal data file
    abnormal_pharmgkb_data = TSV['abnormal_pharmgkb_data']

    ### Check data before load to Mysql
    def check_data(table_name, list_name, format, out_file):
        field_regex = []
        error_num = 0
        for ele in list(SCHEMA[table_name].keys()):
            if '_type' in ele:
                field_regex.append(SCHEMA[table_name][ele])
        for ele in list_name:
            if isinstance(ele, list):
                error_num = error_num + 1
                logger.error('Unexpected data list appears: %s' % list_name)
                # Record the abnormal data into a file
                with open(abnormal_pharmgkb_data, 'a') as abnormal_file:
                    print(list_name, file = abnormal_file)
        for i in range(len(field_regex)):
            matchObj = re.match(field_regex[i], str(list_name[i]))
            if matchObj is None:
                error_num = error_num + 1
                logger.error('Format error occurs: %i %s' % (i, list_name))
                # Record the abnormal data into a file
                with open(abnormal_pharmgkb_data, 'a') as abnormal_file:
                    print(list_name, file = abnormal_file)
        if error_num == 0:
            print(format.format(list_name), file = out_file)
        return(error_num)

    study_type = []
    ratio_stat_type = []
    charact_type = []
    race =[]
    for item in study_parameters:
        temp1 = item[SCHEMA['PharmGKBStudyTypeDic']['Name']]
        study_type.extend(temp1.split(", "))
        ratio_stat_type.append(item[SCHEMA['PharmGKBRatioStatTypeDic']['Name']])
        charact_type.append(item[SCHEMA['PharmGKBCharacteristicsTypeDic']['Name']])
        temp2 = item[SCHEMA['PharmGKBRaceDic']['Name']]
        race.extend(temp2.split(", "))

    for item in clinical_ann_metadata:
        temp3 = item[SCHEMA['PharmGKBRaceDic']['Name']]
        race.extend(temp3.split(", "))

    for i in range(len(race)):
        race[i] = re.sub('in-vitro|in vitro', "In vitro", race[i])

    study_type_dic = list(set(study_type))
    ratio_stat_type_dic = list(set(ratio_stat_type))
    charact_type_dic = list(set(charact_type))
    race_dic = list(set(race))


    # Create PharmGKBStudyTypeDic
    logger.info('Start to create PharmGKBStudyTypeDic')
    tsv1 = TSV['PharmGKBStudyTypeDic']
    total_errors = 0
    with open(tsv1, 'w') as f1:
        for Name in study_type_dic:
            error_num = check_data('PharmGKBStudyTypeDic', [Name], '{0[0]}', f1)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBStudyTypeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv1)
    logger.info('PharmGKBStudyTypeDic records: %i' % row_count)


    # Create PharmGKBRatioStatTypeDic
    logger.info('Start to create PharmGKBRatioStatTypeDic')
    tsv2 = TSV['PharmGKBRatioStatTypeDic']
    total_errors = 0
    with open(tsv2, 'w') as f2:
        for Name in ratio_stat_type_dic:
            error_num = check_data('PharmGKBRatioStatTypeDic', [Name], '{0[0]}', f2)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBRatioStatTypeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv2)
    logger.info('PharmGKBRatioStatTypeDic records: %i' % row_count)


    # Create PharmGKBCharacteristicsTypeDic
    logger.info('Start to create PharmGKBCharacteristicsTypeDic')
    tsv3 = TSV['PharmGKBCharacteristicsTypeDic']
    total_errors = 0
    with open(tsv3, 'w') as f3:
        for Name in charact_type_dic:
            error_num = check_data('PharmGKBCharacteristicsTypeDic', [Name], '{0[0]}', f3)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(3)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBCharacteristicsTypeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv3)
    logger.info('PharmGKBCharacteristicsTypeDic records: %i' % row_count)


    # Create PharmGKBRaceDic
    logger.info('Start to create PharmGKBRaceDic')
    tsv4 = TSV['PharmGKBRaceDic']
    total_errors = 0
    with open(tsv4, 'w') as f4:
        for Name in race_dic:
            error_num = check_data('PharmGKBRaceDic', [Name], '{0[0]}', f4)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(4)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBRaceDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv4)
    logger.info('PharmGKBRaceDic records: %i' % row_count)


    # Create PharmGKBStudyParameter
    logger.info('Start to create PharmGKBStudyParameter')
    tsv5 = TSV['PharmGKBStudyParameter']
    total_errors = 0
    with open(tsv5, 'w') as f5:
        for item in study_parameters:
            StudyParametersID = item[SCHEMA["PharmGKBStudyParameter"]["StudyParametersID"]]
            CaseNum = item[SCHEMA["PharmGKBStudyParameter"]["CaseNum"]]
            ControlNum = item[SCHEMA["PharmGKBStudyParameter"]["ControlNum"]]
            Characteristics = item[SCHEMA["PharmGKBStudyParameter"]["Characteristics"]]
            CharactType = item[SCHEMA['PharmGKBCharacteristicsTypeDic']['Name']]
            CharactTypeID = pymysql_cursor("SELECT ID FROM PharmGKBCharacteristicsTypeDic WHERE Name = '%s';" % CharactType)
            AlleleFreqInCase = item[SCHEMA["PharmGKBStudyParameter"]["AlleleFreqInCase"]]
            AlleleInCase = item[SCHEMA["PharmGKBStudyParameter"]["AlleleInCase"]]
            AlleleFreqInControl = item[SCHEMA["PharmGKBStudyParameter"]["AlleleFreqInControl"]]
            AlleleInControl = item[SCHEMA["PharmGKBStudyParameter"]["AlleleInControl"]]
            PValueOperator = item[SCHEMA["PharmGKBStudyParameter"]["PValueOperator"]]
            PValue = item[SCHEMA["PharmGKBStudyParameter"]["PValue"]]
            RatioStatType = item[SCHEMA['PharmGKBRatioStatTypeDic']['Name']]
            RatioStatTypeID = pymysql_cursor("SELECT ID FROM PharmGKBRatioStatTypeDic WHERE Name = '%s';" % RatioStatType)
            RatioStat = item[SCHEMA["PharmGKBStudyParameter"]['RatioStat']]
            ConfiIntStart = item[SCHEMA["PharmGKBStudyParameter"]['ConfiIntStart']]
            ConfiIntStop = item[SCHEMA["PharmGKBStudyParameter"]['ConfiIntStop']]
            # Check data
            wait_confirm = [StudyParametersID, CaseNum, ControlNum, Characteristics, CharactTypeID, AlleleFreqInCase, AlleleInCase, AlleleFreqInControl, AlleleInControl, PValueOperator, PValue, RatioStatTypeID, RatioStat, ConfiIntStart, ConfiIntStop]
            error_num = check_data('PharmGKBStudyParameter', wait_confirm, '{0[0]}\t{0[1]}\t{0[2]}\t{0[3]}\t{0[4]}\t{0[5]}\t{0[6]}\t{0[7]}\t{0[8]}\t{0[9]}\t{0[10]}\t{0[11]}\t{0[12]}\t{0[13]}\t{0[14]}', f5)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBStudyParameter FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (StudyParametersID, @CaseNum, @ControlNum, Characteristics, CharactTypeID, @AlleleFreqInCase, AlleleInCase, @AlleleFreqInControl, AlleleInControl, PValueOperator, @PValue, RatioStatTypeID, @RatioStat, @ConfiIntStart, @ConfiIntStop) \
    SET CaseNum = CASE WHEN @CaseNum NOT IN ('None') THEN @CaseNum END, \
        ControlNum = CASE WHEN @ControlNum NOT IN ('None') THEN @ControlNum END, \
        AlleleFreqInCase = CASE WHEN @AlleleFreqInCase NOT IN ('None') THEN @AlleleFreqInCase END, \
        AlleleFreqInControl = CASE WHEN @AlleleFreqInControl NOT IN ('None') THEN @AlleleFreqInControl END, \
        PValue = CASE WHEN @PValue NOT IN ('None') THEN @PValue END, \
        RatioStat = CASE WHEN @RatioStat NOT IN ('None') THEN @RatioStat END, \
        ConfiIntStart = CASE WHEN @ConfiIntStart NOT IN ('None') THEN @ConfiIntStart END, \
        ConfiIntStop = CASE WHEN @ConfiIntStop NOT IN ('None') THEN @ConfiIntStop END;" % tsv5)
    logger.info('PharmGKBStudyParameter records: %i' % row_count)


    # Create PharmGKBStuParaHasRace
    logger.info('Start to create PharmGKBStuParaHasRace')
    tsv6 = TSV['PharmGKBStuParaHasRace']
    total_errors = 0
    with open(tsv6, 'w') as f6:
        for item in study_parameters:
            StudyParametersID = item[SCHEMA["PharmGKBStudyParameter"]["StudyParametersID"]]
            StuParaID = pymysql_cursor("SELECT ID FROM PharmGKBStudyParameter WHERE StudyParametersID = '%s';" % StudyParametersID)
            races = item[SCHEMA['PharmGKBRaceDic']['Name']]
            races = re.sub('in-vitro|in vitro', "In vitro", races)
            for Race in races.split(', '):
                RaceID = pymysql_cursor("SELECT ID FROM PharmGKBRaceDic WHERE Name = '%s';" % Race)
                # Check data
                wait_confirm = [StuParaID, RaceID]
                error_num = check_data('PharmGKBStuParaHasRace', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f6)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBStuParaHasRace FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (StuParaID, RaceID);" % tsv6)
    logger.info('PharmGKBStuParaHasRace records: %i' % row_count)


    # Create PharmGKBStuParaHasStuType
    logger.info('Start to create PharmGKBStuParaHasStuType')
    tsv7 = TSV['PharmGKBStuParaHasStuType']
    total_errors = 0
    with open(tsv7, 'w') as f7:
        for item in study_parameters:
            StudyParametersID = item[SCHEMA["PharmGKBStudyParameter"]["StudyParametersID"]]
            StuParaID = pymysql_cursor("SELECT ID FROM PharmGKBStudyParameter WHERE StudyParametersID = '%s';" % StudyParametersID)
            study_types = item[SCHEMA['PharmGKBStudyTypeDic']['Name']]
            for StuType in study_types.split(', '):
                StuTypeID = pymysql_cursor("SELECT ID FROM PharmGKBStudyTypeDic WHERE Name = '%s';" % StuType)
                # Check data
                wait_confirm = [StuParaID, StuTypeID]
                error_num = check_data('PharmGKBStuParaHasStuType', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f7)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBStuParaHasStuType FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (StuParaID, StuTypeID);" % tsv7)
    logger.info('PharmGKBStuParaHasStuType records: %i' % row_count)


    significance = []
    chromosome = []
    clinical_ann_type = []
    for item in var_ann:
        significance.append(item[SCHEMA["PharmGKBSignificanceDic"]["Name"]])
        chromosome.append(item[SCHEMA["PharmGKBChromosomeDic"]["Name"]])
        temp = re.split(',|;', item[SCHEMA["PharmGKBClinAnnTypeDic"]["Name"][0]].strip())
        for ele in temp:
            typ = ele.strip('"')
            if typ != "":
                typ = typ[0].upper() + typ[1:]
                clinical_ann_type.append(typ)
            else:
                clinical_ann_type.append(typ)

    for item in clinical_ann_metadata:
        temp = re.split(',|;', item[SCHEMA["PharmGKBClinAnnTypeDic"]["Name"][1]].strip())
        for ele in temp:
            typ = ele.strip('"')
            if typ != "":
                typ = typ[0].upper() + typ[1:]
                clinical_ann_type.append(typ)
            else:
                clinical_ann_type.append(typ)

    significance_dic = list(set(significance))
    chromosome_dic = list(set(chromosome))
    clinical_ann_type_dic = list(set(clinical_ann_type))


    # Create PharmGKBSignificanceDic
    logger.info('Start to create PharmGKBSignificanceDic')
    tsv8 = TSV['PharmGKBSignificanceDic']
    total_errors = 0
    with open(tsv8, 'w') as f8:
        for Name in significance_dic:
            error_num = check_data('PharmGKBSignificanceDic', [Name], '{0[0]}', f8)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBSignificanceDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv8)
    logger.info('PharmGKBSignificanceDic records: %i' % row_count)


    # Create PharmGKBChromosomeDic
    logger.info('Start to create PharmGKBChromosomeDic')
    tsv9 = TSV['PharmGKBChromosomeDic']
    total_errors = 0
    with open(tsv9, 'w') as f9:
        for Name in chromosome_dic:
            error_num = check_data('PharmGKBChromosomeDic', [Name], '{0[0]}', f9)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBChromosomeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv9)
    logger.info('PharmGKBChromosomeDic records: %i' % row_count)


    # Create PharmGKBClinAnnTypeDic
    logger.info('Start to create PharmGKBClinAnnTypeDic')
    tsv10 = TSV['PharmGKBClinAnnTypeDic']
    total_errors = 0
    with open(tsv10, 'w') as f10:
        for Name in clinical_ann_type_dic:
            error_num = check_data('PharmGKBClinAnnTypeDic', [Name], '{0[0]}', f10)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnTypeDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv10)
    logger.info('PharmGKBClinAnnTypeDic records: %i' % row_count)


    # Create PharmGKBVarAnn
    logger.info('Start to create PharmGKBVarAnn')
    tsv11 = TSV['PharmGKBVarAnn']
    total_errors = 0
    with open(tsv11, 'w') as f11:
        for item in var_ann:
            AnnotationID = item[SCHEMA["PharmGKBVarAnn"]["AnnotationID"]]
            Variant = item[SCHEMA["PharmGKBVarAnn"]["Variant"]]
            Gene = re.sub('"', '', item[SCHEMA["PharmGKBVarAnn"]["Gene"]])
            Drug = re.sub('"', '', item[SCHEMA["PharmGKBVarAnn"]["Drug"]])
            PMID = item[SCHEMA["PharmGKBVarAnn"]["PMID"]]
            Significance = item[SCHEMA["PharmGKBVarAnn"]["SignificanceID"]]
            SignificanceID = pymysql_cursor("SELECT ID FROM PharmGKBSignificanceDic WHERE Name = '%s';" % Significance)
            Notes = item[SCHEMA["PharmGKBVarAnn"]["Notes"]]
            Sentence = item[SCHEMA["PharmGKBVarAnn"]["Sentence"]]
            Allele = item[SCHEMA["PharmGKBVarAnn"]["Allele"]]
            Chromosome = item[SCHEMA["PharmGKBVarAnn"]["ChrID"]]
            ChrID = pymysql_cursor("SELECT ID FROM PharmGKBChromosomeDic WHERE Name = '%s';" % Chromosome)
            # Check data
            wait_confirm = [AnnotationID, Variant, Gene, Drug, PMID, SignificanceID, Notes, Sentence, Allele, ChrID]
            error_num = check_data('PharmGKBVarAnn', wait_confirm, '{0[0]:.0f}\t{0[1]}\t{0[2]}\t{0[3]}\t{0[4]}\t{0[5]:.0f}\t{0[6]}\t{0[7]}\t{0[8]}\t{0[9]:.0f}', f11)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor('LOAD DATA LOCAL INFILE "%s" INTO TABLE PharmGKBVarAnn FIELDS TERMINATED BY "\\t" LINES TERMINATED BY "\\n" (AnnotationID, Variant, Gene, Drug, PMID, SignificanceID, Notes, Sentence, Allele, ChrID);' % tsv11)
    logger.info('PharmGKBVarAnn records: %i' % row_count)


    # Create PharmGKBVarAnnHasClinAnnType
    logger.info('Start to create PharmGKBVarAnnHasClinAnnType')
    tsv12 = TSV['PharmGKBVarAnnHasClinAnnType']
    total_errors = 0
    with open(tsv12, 'w') as f12:
        for item in var_ann:
            AnnotationID = item[SCHEMA["PharmGKBVarAnn"]["AnnotationID"]]
            VarAnnID = pymysql_cursor("SELECT ID FROM PharmGKBVarAnn WHERE AnnotationID = '%s';" % AnnotationID)
            ClinAnnType_temp = re.split(',|;', item[SCHEMA["PharmGKBClinAnnTypeDic"]["Name"][0]].strip())
            for ele in ClinAnnType_temp:
                typ = ele.strip('"')
                if typ != "":
                    ClinAnnType = typ[0].upper() + typ[1:]
                else:
                    ClinAnnType = typ
                ClinAnnTypeID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnnTypeDic WHERE Name = '%s';" % ClinAnnType)
                # Check data
                wait_confirm = [VarAnnID, ClinAnnTypeID]
                error_num = check_data('PharmGKBVarAnnHasClinAnnType', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f12)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBVarAnnHasClinAnnType FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (VarAnnID, ClinAnnTypeID);" % tsv12)
    logger.info('PharmGKBVarAnnHasClinAnnType records: %i' % row_count)


    # Create PharmGKBVarAnnHasStuPara
    logger.info('Start to create PharmGKBVarAnnHasStuPara')
    tsv13 = TSV['PharmGKBVarAnnHasStuPara']
    total_errors = 0
    with open(tsv13, 'w') as f13:
        for item in var_ann:
            AnnotationID = item[SCHEMA["PharmGKBVarAnn"]["AnnotationID"]]
            VarAnnID = pymysql_cursor("SELECT ID FROM PharmGKBVarAnn WHERE AnnotationID = '%s';" % AnnotationID)
            StuParaIDs = re.sub('"|\'', '', item[SCHEMA["PharmGKBVarAnnHasStuPara"]["StuParaID"]])
            if StuParaIDs != '':
                for stu_para_id in StuParaIDs.split(','):
                    StuParaID = pymysql_cursor("SELECT ID FROM PharmGKBStudyParameter WHERE StudyParametersID = '%s';" % stu_para_id.strip())
                    # Check data
                    wait_confirm = [VarAnnID, StuParaID]
                    error_num = check_data('PharmGKBVarAnnHasStuPara', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f13)
                    total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBVarAnnHasStuPara FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (VarAnnID, StuParaID);" % tsv13)
    logger.info('PharmGKBVarAnnHasStuPara records: %i' % row_count)


    # Create PharmGKBLevelOfEvidenceDic
    evidence_level = []
    for item in clinical_ann_metadata:
        evidence_level.append(item[SCHEMA["PharmGKBLevelOfEvidenceDic"]["Name"]])

    evidence_level_dic = list(set(evidence_level))
    logger.info('Start to create PharmGKBLevelOfEvidenceDic')
    tsv14 = TSV['PharmGKBLevelOfEvidenceDic']
    total_errors = 0
    with open(tsv14, 'w') as f14:
        for Name in evidence_level_dic:
            error_num = check_data('PharmGKBLevelOfEvidenceDic', [Name], '{0[0]}', f14)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBLevelOfEvidenceDic FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv14)
    logger.info('PharmGKBLevelOfEvidenceDic records: %i' % row_count)


    # Create PharmGKBPhenoGeno
    logger.info('Start to create PharmGKBPhenoGeno')
    tsv15 = TSV['PharmGKBPhenoGeno']
    total_errors = 0
    with open(tsv15, 'w') as f15:
        for item in clinical_ann:
            PhenoGenotypeID = item[SCHEMA["PharmGKBPhenoGeno"]["PhenoGenotypeID"]]
            Genotype = item[SCHEMA["PharmGKBPhenoGeno"]["Genotype"]]
            ClinicalPhenotype = item[SCHEMA["PharmGKBPhenoGeno"]["ClinicalPhenotype"]]
            # Check data
            wait_confirm = [PhenoGenotypeID, Genotype, ClinicalPhenotype]
            error_num = check_data('PharmGKBPhenoGeno', wait_confirm, '{0[0]:.0f}\t{0[1]}\t{0[2]}', f15)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBPhenoGeno FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (PhenoGenotypeID, Genotype, ClinicalPhenotype);" % tsv15)
    logger.info('PharmGKBPhenoGeno records: %i' % row_count)


    all_genes = []
    all_vars = []
    all_drugs = []
    all_diseases = []
    for item in clinical_ann_metadata:
        Gene = item[SCHEMA["PharmGKBGene"]["Name"]]
        Variant = item[SCHEMA["PharmGKBVariant"]["Name"]]
        Drug = item[SCHEMA["PharmGKBDrug"]["Name"]]
        Disease = item[SCHEMA["PharmGKBDisease"]["Name"]]
        all_genes.extend(Gene.split(";"))
        if Variant == "G6PD B (wildtype), G6PD Canton, Taiwan-Hakka, Gifu-like, Agrigento-like":
            Variant = "G6PD B (wildtype), G6PD Canton, G6PD Taiwan-Hakka, G6PD Gifu-like, G6PD Agrigento-like"
        if Variant == "G6PD B (wildtype), G6PD Mediterranean, Dallas, Panama' Sassari, Cagliari, Birmingham":
            Variant = "G6PD B (wildtype), G6PD Mediterranean, G6PD Dallas, G6PD Panama' Sassari, G6PD Cagliari, G6PD Birmingham"
        all_vars.extend(Variant.split(", "))
        all_drugs.extend(Drug.split(";"))
        all_diseases.extend(Disease.split(";"))

    all_genes = list(set(all_genes))
    all_vars = list(set(all_vars))
    all_drugs = list(set(all_drugs))
    all_diseases = list(set(all_diseases))

    # Create PharmGKBGene
    logger.info('Start to create PharmGKBGene')
    tsv21 = TSV['PharmGKBGene']
    total_errors = 0
    with open(tsv21, 'w') as f21:
        for item in all_genes:
            paid = re.findall(r'[(PA](\d*?)[)]', item)
            if len(paid) > 1:
                print(item)
            if paid:
                PAID = 'PA' + paid[0]
            else:
                PAID = ''
            to_replace = '(' + PAID + ')'
            Name = item.replace(to_replace, '').strip()
            # Check data
            wait_confirm = [Name, PAID]
            error_num = check_data('PharmGKBGene', wait_confirm, '{0[0]}\t{0[1]}', f21)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBGene FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name, PAID);" % tsv21)
    logger.info('PharmGKBGene records: %i' % row_count)

    # Create PharmGKBVariant
    logger.info('Start to create PharmGKBVariant')
    tsv22 = TSV['PharmGKBVariant']
    total_errors = 0
    with open(tsv22, 'w') as f22:
        for item in all_vars:
            Name = item.strip()
            # Check data
            wait_confirm = [Name]
            error_num = check_data('PharmGKBVariant', wait_confirm, '{0[0]}', f22)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBVariant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name);" % tsv22)
    logger.info('PharmGKBVariant records: %i' % row_count)

    # Create PharmGKBDrug
    logger.info('Start to create PharmGKBDrug')
    tsv23 = TSV['PharmGKBDrug']
    total_errors = 0
    with open(tsv23, 'w') as f23:
        for item in all_drugs:
            paid = re.findall(r'[(PA](\d*?)[)]', item)
            if len(paid) > 1:
                print(item)
            if paid:
                PAID = 'PA' + paid[0]
            else:
                PAID = ''
            to_replace = '(' + PAID + ')'
            Name = item.replace(to_replace, '').strip()
            # Check data
            wait_confirm = [Name, PAID]
            error_num = check_data('PharmGKBDrug', wait_confirm, '{0[0]}\t{0[1]}', f23)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBDrug FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name, PAID);" % tsv23)
    logger.info('PharmGKBDrug records: %i' % row_count)

    # Create PharmGKBDisease
    logger.info('Start to create PharmGKBDisease')
    tsv24 = TSV['PharmGKBDisease']
    total_errors = 0
    with open(tsv24, 'w') as f24:
        for item in all_diseases:
            paid = re.findall(r'[(PA](\d*?)[)]', item)
            if len(paid) > 1:
                print(item)
            if paid:
                PAID = 'PA' + paid[0]
            else:
                PAID = ''
            to_replace = '(' + PAID + ')'
            Name = item.replace(to_replace, '').strip()
            # Check data
            wait_confirm = [Name, PAID]
            error_num = check_data('PharmGKBDisease', wait_confirm, '{0[0]}\t{0[1]}', f24)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBDisease FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (Name, PAID);" % tsv24)
    logger.info('PharmGKBDisease records: %i' % row_count)


    # Create PharmGKBClinAnn
    logger.info('Start to create PharmGKBClinAnn')
    tsv16 = TSV['PharmGKBClinAnn']
    total_errors = 0
    with open(tsv16, 'w') as f16:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            EvidenceLevel = item[SCHEMA["PharmGKBLevelOfEvidenceDic"]["Name"]]
            EvidenceLevelID = pymysql_cursor("SELECT ID FROM PharmGKBLevelOfEvidenceDic WHERE Name = '%s';" % EvidenceLevel)
            EvidenceCount = item[SCHEMA["PharmGKBClinAnn"]["EvidenceCount"]]
            # Check data
            wait_confirm = [ClinicalAnnotationID, EvidenceLevelID, EvidenceCount]
            error_num = check_data('PharmGKBClinAnn', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}\t{0[2]:.0f}', f16)
            total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnn FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinicalAnnotationID, EvidenceLevelID, EvidenceCount);" % tsv16)
    logger.info('PharmGKBClinAnn records: %i' % row_count)


    # Create PharmGKBClinAnnHasGene
    logger.info('Start to create PharmGKBClinAnnHasGene')
    tsv24 = TSV['PharmGKBClinAnnHasGene']
    total_errors = 0
    with open(tsv24, 'w') as f24:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            Genes = item[SCHEMA["PharmGKBGene"]["Name"]].split(";")
            for item in Genes:
                paid = re.findall(r'[(PA](\d*?)[)]', item)
                if paid:
                    PAID = 'PA' + paid[0]
                else:
                    PAID = ''
                GeneID = pymysql_cursor("SELECT ID FROM PharmGKBGene WHERE PAID = '%s';" % PAID)
                # Check data
                wait_confirm = [ClinAnnID, GeneID]
                error_num = check_data('PharmGKBClinAnnHasGene', wait_confirm, '{0[0]}\t{0[1]}', f24)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasGene FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, GeneID);" % tsv24)
    logger.info('PharmGKBClinAnnHasGene records: %i' % row_count)


    # Create PharmGKBClinAnnHasVariant
    logger.info('Start to create PharmGKBClinAnnHasVariant')
    tsv24 = TSV['PharmGKBClinAnnHasVariant']
    total_errors = 0
    with open(tsv24, 'w') as f24:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            Variants = item[SCHEMA["PharmGKBVariant"]["Name"]]
            if Variants == "G6PD B (wildtype), G6PD Canton, Taiwan-Hakka, Gifu-like, Agrigento-like":
                Variants = "G6PD B (wildtype), G6PD Canton, G6PD Taiwan-Hakka, G6PD Gifu-like, G6PD Agrigento-like"
            elif Variants == "G6PD B (wildtype), G6PD Mediterranean, Dallas, Panama' Sassari, Cagliari, Birmingham":
                Variants = "G6PD B (wildtype), G6PD Mediterranean, G6PD Dallas, G6PD Panama' Sassari, G6PD Cagliari, G6PD Birmingham"
            Variants = Variants.split(", ")
            for item in Variants:
                VariantID = pymysql_cursor('SELECT ID FROM PharmGKBVariant WHERE Name = "%s";' % item.strip())
                # Check data
                wait_confirm = [ClinAnnID, VariantID]
                error_num = check_data('PharmGKBClinAnnHasVariant', wait_confirm, '{0[0]}\t{0[1]}', f24)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasVariant FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, VariantID);" % tsv24)
    logger.info('PharmGKBClinAnnHasVariant records: %i' % row_count)


    # Create PharmGKBClinAnnHasDrug
    logger.info('Start to create PharmGKBClinAnnHasDrug')
    tsv24 = TSV['PharmGKBClinAnnHasDrug']
    total_errors = 0
    with open(tsv24, 'w') as f24:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            Drugs = item[SCHEMA["PharmGKBDrug"]["Name"]].split(";")
            for item in Drugs:
                paid = re.findall(r'[(PA](\d*?)[)]', item)
                if paid:
                    PAID = 'PA' + paid[0]
                else:
                    PAID = ''
                DrugID = pymysql_cursor("SELECT ID FROM PharmGKBDrug WHERE PAID = '%s';" % PAID)
                # Check data
                wait_confirm = [ClinAnnID, DrugID]
                error_num = check_data('PharmGKBClinAnnHasDrug', wait_confirm, '{0[0]}\t{0[1]}', f24)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasDrug FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, DrugID);" % tsv24)
    logger.info('PharmGKBClinAnnHasDrug records: %i' % row_count)


    # Create PharmGKBClinAnnHasDisease
    logger.info('Start to create PharmGKBClinAnnHasDisease')
    tsv24 = TSV['PharmGKBClinAnnHasDisease']
    total_errors = 0
    with open(tsv24, 'w') as f24:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            Diseases = item[SCHEMA["PharmGKBDisease"]["Name"]].split(";")
            for item in Diseases:
                paid = re.findall(r'[(PA](\d*?)[)]', item)
                if paid:
                    PAID = 'PA' + paid[0]
                else:
                    PAID = ''
                DiseaseID = pymysql_cursor("SELECT ID FROM PharmGKBDisease WHERE PAID = '%s';" % PAID)
                # Check data
                wait_confirm = [ClinAnnID, DiseaseID]
                error_num = check_data('PharmGKBClinAnnHasDisease', wait_confirm, '{0[0]}\t{0[1]}', f24)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasDisease FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, DiseaseID);" % tsv24)
    logger.info('PharmGKBClinAnnHasDisease records: %i' % row_count)


    # Create PharmGKBClinAnnHasPhenoGeno
    logger.info('Start to create PharmGKBClinAnnHasPhenoGeno')
    tsv17 = TSV['PharmGKBClinAnnHasPhenoGeno']
    total_errors = 0
    with open(tsv17, 'w') as f17:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            PhenoGenoIDs = item[SCHEMA["PharmGKBClinAnnHasPhenoGeno"]["PhenoGenoID"]]
            for PhenoGenoID in PhenoGenoIDs.split(';'):
                PhenoGenoID = pymysql_cursor("SELECT ID FROM PharmGKBPhenoGeno WHERE PhenoGenotypeID = '%s';" % PhenoGenoID)
                # Check data
                wait_confirm = [ClinAnnID, PhenoGenoID]
                error_num = check_data('PharmGKBClinAnnHasPhenoGeno', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f17)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasPhenoGeno FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, PhenoGenoID);" % tsv17)
    logger.info('PharmGKBClinAnnHasPhenoGeno records: %i' % row_count)


    # Create PharmGKBClinAnnHasRace
    logger.info('Start to create PharmGKBClinAnnHasRace')
    tsv18 = TSV['PharmGKBClinAnnHasRace']
    total_errors = 0
    with open(tsv18, 'w') as f18:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            races = item[SCHEMA['PharmGKBRaceDic']['Name']]
            races = re.sub('in-vitro|in vitro', "In vitro", races)
            for Race in races.split(', '):
                RaceID = pymysql_cursor("SELECT ID FROM PharmGKBRaceDic WHERE Name = '%s';" % Race)
                # Check data
                wait_confirm = [ClinAnnID, RaceID]
                error_num = check_data('PharmGKBClinAnnHasRace', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f18)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasRace FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, RaceID);" % tsv18)
    logger.info('PharmGKBClinAnnHasRace records: %i' % row_count)


    # Create PharmGKBClinAnnHasClinAnnType
    logger.info('Start to create PharmGKBClinAnnHasClinAnnType')
    tsv19 = TSV['PharmGKBClinAnnHasClinAnnType']
    total_errors = 0
    with open(tsv19, 'w') as f19:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            ClinAnnType_temp = re.split(',|;', item[SCHEMA["PharmGKBClinAnnTypeDic"]["Name"][1]].strip())
            for ele in ClinAnnType_temp:
                typ = ele.strip('"')
                if typ != "":
                    ClinAnnType = typ[0].upper() + typ[1:]
                else:
                    ClinAnnType = typ
                ClinAnnTypeID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnnTypeDic WHERE Name = '%s';" % ClinAnnType)
                # Check data
                wait_confirm = [ClinAnnID, ClinAnnTypeID]
                error_num = check_data('PharmGKBClinAnnHasClinAnnType', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f19)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasClinAnnType FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, ClinAnnTypeID);" % tsv19)
    logger.info('PharmGKBClinAnnHasClinAnnType records: %i' % row_count)


    # Create PharmGKBClinAnnHasVarAnn
    logger.info('Start to create PharmGKBClinAnnHasVarAnn')
    tsv20 = TSV['PharmGKBClinAnnHasVarAnn']
    total_errors = 0
    with open(tsv20, 'w') as f20:
        for item in clinical_ann_metadata:
            ClinicalAnnotationID = item[SCHEMA["PharmGKBClinAnn"]["ClinicalAnnotationID"]]
            ClinAnnID = pymysql_cursor("SELECT ID FROM PharmGKBClinAnn WHERE ClinicalAnnotationID = '%s';" % ClinicalAnnotationID)
            AnnotationIDs = item[SCHEMA["PharmGKBClinAnnHasVarAnn"]["VarAnnID"]]
            for var_ann_id in AnnotationIDs.split(';'):
                VarAnnID = pymysql_cursor("SELECT ID FROM PharmGKBVarAnn WHERE AnnotationID = '%s';" % var_ann_id)
                # Check data
                wait_confirm = [ClinAnnID, VarAnnID]
                error_num = check_data('PharmGKBClinAnnHasVarAnn', wait_confirm, '{0[0]:.0f}\t{0[1]:.0f}', f20)
                total_errors = total_errors + error_num
        if total_errors != 0:
            print("error")
            #sys.exit(2)

    # Load data
    row_count = pymysql_cursor("LOAD DATA LOCAL INFILE '%s' INTO TABLE PharmGKBClinAnnHasVarAnn FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (ClinAnnID, VarAnnID);" % tsv20)
    logger.info('PharmGKBClinAnnHasVarAnn records: %i' % row_count)


### Logging
def load_logging_cfg():
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

if __name__ == '__main__':
    main(sys.argv[1:])
    exit(0)