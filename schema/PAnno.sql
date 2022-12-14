-- MySQL Script generated by MySQL Workbench
-- Wed Nov 23 10:40:47 2022
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';


-- -----------------------------------------------------
-- Schema PAnno
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `PAnno` DEFAULT CHARACTER SET utf8 ;
USE `PAnno` ;

-- -----------------------------------------------------
-- Table `PAnno`.`CharacteristicsTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`CharacteristicsTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(32) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB
COMMENT = 'Charicteristic type dictionary';


-- -----------------------------------------------------
-- Table `PAnno`.`Gene`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`Gene` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Symbol` VARCHAR(32) NOT NULL,
  `Name` VARCHAR(256) NULL DEFAULT NULL,
  `PAID` VARCHAR(16) NULL DEFAULT NULL,
  `EntrezID` VARCHAR(32) NULL DEFAULT NULL,
  `HGNCID` VARCHAR(32) NULL DEFAULT NULL,
  `EnsemblID` VARCHAR(512) NULL DEFAULT NULL,
  `Chr` VARCHAR(8) NULL DEFAULT NULL,
  `GRCh37Start` VARCHAR(16) NULL DEFAULT NULL,
  `GRCh37End` VARCHAR(16) NULL DEFAULT NULL,
  `GRCh38Start` VARCHAR(16) NULL DEFAULT NULL,
  `GRCh38End` VARCHAR(16) NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `PAID_UNIQUE` (`PAID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`Variant`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`Variant` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL COMMENT 'Variant(rsID) or haplotype',
  `NCID` VARCHAR(32) NULL DEFAULT NULL,
  `PAID` VARCHAR(16) NULL DEFAULT NULL,
  `Type` VARCHAR(16) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`Drug`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`Drug` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(128) NOT NULL,
  `PAID` VARCHAR(16) NOT NULL,
  `GenericName` VARCHAR(1024) NULL DEFAULT NULL,
  `TradeName` VARCHAR(128) NULL DEFAULT NULL,
  `CrossRef` TEXT NULL DEFAULT NULL,
  `Type` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `PAID_UNIQUE` (`PAID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`PhenotypeCategoryDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`PhenotypeCategoryDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`Phenotype`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`Phenotype` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(128) NOT NULL,
  `PAID` VARCHAR(16) NOT NULL,
  `Type` VARCHAR(32) NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `PAID_UNIQUE` (`PAID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`ClinAnn`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`ClinAnn` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'ClinicalAnn',
  `CAID` INT(16) UNSIGNED NOT NULL COMMENT 'Clinical Annotation ID',
  `Gene` VARCHAR(64) NOT NULL,
  `Variant` VARCHAR(512) NOT NULL,
  `Allele1` VARCHAR(64) NOT NULL,
  `Allele2` VARCHAR(64) NULL DEFAULT NULL,
  `Annotation1` TEXT NOT NULL,
  `Annotation2` TEXT NULL DEFAULT NULL,
  `Function1` VARCHAR(64) NOT NULL,
  `Function2` VARCHAR(64) NULL DEFAULT NULL,
  `Score1` FLOAT(8) NULL DEFAULT NULL,
  `Score2` FLOAT(8) NULL DEFAULT NULL,
  `CPICPhenotype` VARCHAR(64) NULL DEFAULT NULL,
  `PAnnoPhenotype` VARCHAR(64) NOT NULL,
  `Drug` VARCHAR(256) NOT NULL,
  `Phenotypes` VARCHAR(256) NOT NULL,
  `EvidenceLevel` VARCHAR(4) NOT NULL,
  `LevelOverride` TEXT NULL DEFAULT NULL,
  `LevelModifier` VARCHAR(32) NULL DEFAULT NULL,
  `Score` VARCHAR(32) NOT NULL,
  `PMIDCount` INT(4) UNSIGNED NOT NULL,
  `EvidenceCount` INT(4) NOT NULL,
  `Specialty` VARCHAR(16) NULL DEFAULT NULL,
  `PhenotypeCategory` VARCHAR(16) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB
COMMENT = 'Original record of clinical annotation.';


-- -----------------------------------------------------
-- Table `PAnno`.`ClinAnnHasPhenotype`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`ClinAnnHasPhenotype` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `ClinAnnID` INT(8) UNSIGNED NOT NULL,
  `PhenotypeID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnnHasPhenotype_ClinAnn1_idx` (`ClinAnnID` ASC),
  INDEX `fk_ClinAnnHasPhenotype_Phenotype1_idx` (`PhenotypeID` ASC),
  CONSTRAINT `fk_ClinAnnHasPhenotype_ClinAnn1`
    FOREIGN KEY (`ClinAnnID`)
    REFERENCES `PAnno`.`ClinAnn` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnnHasPhenotype_Phenotype1`
    FOREIGN KEY (`PhenotypeID`)
    REFERENCES `PAnno`.`Phenotype` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`ClinAnnHasVariant`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`ClinAnnHasVariant` (
  `ID` INT(8) NOT NULL AUTO_INCREMENT,
  `ClinAnnID` INT(8) UNSIGNED NOT NULL,
  `VariantID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnnHasVariant_ClinAnn1_idx` (`ClinAnnID` ASC),
  INDEX `fk_ClinAnnHasVariant_Variant1_idx` (`VariantID` ASC),
  CONSTRAINT `fk_ClinAnnHasVariant_ClinAnn1`
    FOREIGN KEY (`ClinAnnID`)
    REFERENCES `PAnno`.`ClinAnn` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnnHasVariant_Variant1`
    FOREIGN KEY (`VariantID`)
    REFERENCES `PAnno`.`Variant` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`ClinAnnHasDrug`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`ClinAnnHasDrug` (
  `ID` INT(8) NOT NULL AUTO_INCREMENT,
  `ClinAnnID` INT(8) UNSIGNED NOT NULL,
  `DrugID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnnHasDrug_ClinAnn1_idx` (`ClinAnnID` ASC),
  INDEX `fk_ClinAnnHasDrug_Drug1_idx` (`DrugID` ASC),
  CONSTRAINT `fk_ClinAnnHasDrug_ClinAnn1`
    FOREIGN KEY (`ClinAnnID`)
    REFERENCES `PAnno`.`ClinAnn` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnnHasDrug_Drug1`
    FOREIGN KEY (`DrugID`)
    REFERENCES `PAnno`.`Drug` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`RaceDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`RaceDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB
COMMENT = 'dictionary: races';


-- -----------------------------------------------------
-- Table `PAnno`.`EvidenceTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`EvidenceTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`VarAnn`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`VarAnn` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `AnnotationID` INT(16) UNSIGNED NOT NULL,
  `VariantOrHaplotype` VARCHAR(512) NOT NULL,
  `Gene` VARCHAR(512) NOT NULL,
  `Drug` VARCHAR(512) NOT NULL,
  `PMID` VARCHAR(64) NULL DEFAULT NULL COMMENT 'PMID of article describing the relation',
  `CategoryID` INT(8) UNSIGNED NOT NULL,
  `Significance` VARCHAR(16) NOT NULL COMMENT 'Foreign key -> PharmGKBSignificance.ID',
  `Note` TEXT NULL DEFAULT NULL COMMENT 'Full describtion of the relation in different heterozygotes/homozygotes',
  `Sentence` TEXT NULL DEFAULT NULL COMMENT 'Summary of Notes',
  `Allele` VARCHAR(256) NULL DEFAULT NULL COMMENT 'Allele studied in the article',
  `Specialty` VARCHAR(64) NULL DEFAULT NULL,
  `TypeID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  FULLTEXT INDEX `Notes` (`Note`),
  FULLTEXT INDEX `Sentence` (`Sentence`),
  INDEX `fk_VarAnn_EvidenceTypeDic1_idx` (`TypeID` ASC),
  INDEX `fk_VarAnn_PhenotypeCategoryDic1_idx` (`CategoryID` ASC),
  CONSTRAINT `fk_VarAnn_EvidenceTypeDic1`
    FOREIGN KEY (`TypeID`)
    REFERENCES `PAnno`.`EvidenceTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_VarAnn_PhenotypeCategoryDic1`
    FOREIGN KEY (`CategoryID`)
    REFERENCES `PAnno`.`PhenotypeCategoryDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'variation annotation main table';


-- -----------------------------------------------------
-- Table `PAnno`.`RatioStatTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`RatioStatTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(8) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB
COMMENT = 'Odds ratio statistic type dictionary';


-- -----------------------------------------------------
-- Table `PAnno`.`StudyParameter`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`StudyParameter` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `StudyParametersID` INT(16) UNSIGNED NOT NULL,
  `CaseNum` INT(8) NULL DEFAULT '0' COMMENT 'Number of patients in case group',
  `ControlNum` INT(8) NULL DEFAULT '0' COMMENT 'Number of patients in control group',
  `Characteristics` TEXT NULL DEFAULT NULL COMMENT 'Study characteristics',
  `CharactTypeID` INT(8) UNSIGNED NOT NULL DEFAULT '0' COMMENT 'Foreign key -> PharmGKBCharType.ID',
  `AlleleFreqInCase` FLOAT NULL DEFAULT '0' COMMENT 'Allele frequency in case group',
  `AlleleInCase` VARCHAR(128) NULL DEFAULT NULL COMMENT 'Alleles in case group',
  `AlleleFreqInControl` FLOAT NULL DEFAULT NULL COMMENT 'Allele frequency in control group',
  `AlleleInControl` VARCHAR(128) NULL DEFAULT NULL COMMENT 'Alleles in control group',
  `PValueOperator` VARCHAR(64) NULL DEFAULT NULL COMMENT '= / > /< in natural language',
  `PValue` FLOAT NULL DEFAULT NULL COMMENT 'p-value',
  `RatioStatTypeID` INT(8) UNSIGNED NOT NULL DEFAULT '0' COMMENT 'Foreign key -> PharmGKBRatioStatType.ID',
  `RatioStat` FLOAT NULL DEFAULT NULL COMMENT 'Odds ratio result',
  `ConfiIntStart` FLOAT NULL DEFAULT NULL COMMENT 'Lower cutoff of confidence interval',
  `ConfiIntStop` FLOAT NULL DEFAULT NULL COMMENT 'Upper cufoff of confidence interval',
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_PharmGKBStudyParameter_PharmGKBCharType1_idx` (`CharactTypeID` ASC),
  INDEX `fk_PharmGKBStudyParameter_PharmGKBRatioStatType1_idx` (`RatioStatTypeID` ASC),
  CONSTRAINT `fk_PharmGKBStudyParameter_PharmGKBCharType1`
    FOREIGN KEY (`CharactTypeID`)
    REFERENCES `PAnno`.`CharacteristicsTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBStudyParameter_PharmGKBRatioStatType1`
    FOREIGN KEY (`RatioStatTypeID`)
    REFERENCES `PAnno`.`RatioStatTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'Study parameters';


-- -----------------------------------------------------
-- Table `PAnno`.`StuParaHasRace`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`StuParaHasRace` (
  `StuParaID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `RaceID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBRace1_idx` (`RaceID` ASC),
  INDEX `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBStudyPar_idx` (`StuParaID` ASC),
  CONSTRAINT `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBRace1`
    FOREIGN KEY (`RaceID`)
    REFERENCES `PAnno`.`RaceDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBStudyParam1`
    FOREIGN KEY (`StuParaID`)
    REFERENCES `PAnno`.`StudyParameter` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'PharmGKBStudyParameter to PharmGkBRace';


-- -----------------------------------------------------
-- Table `PAnno`.`StudyTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`StudyTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB
COMMENT = 'Study type dictionary';


-- -----------------------------------------------------
-- Table `PAnno`.`StuParaHasStuType`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`StuParaHasStuType` (
  `StuTypeID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `StuParaID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStu_idx` (`StuParaID` ASC),
  INDEX `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStu_idx1` (`StuTypeID` ASC),
  CONSTRAINT `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStudy1`
    FOREIGN KEY (`StuTypeID`)
    REFERENCES `PAnno`.`StudyTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStudy2`
    FOREIGN KEY (`StuParaID`)
    REFERENCES `PAnno`.`StudyParameter` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'PharmGKBStudyType to PharmGKBStudyParamet';


-- -----------------------------------------------------
-- Table `PAnno`.`VarAnnHasStuPara`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`VarAnnHasStuPara` (
  `VarAnnID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `StuParaID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBStudyP_idx` (`StuParaID` ASC),
  INDEX `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBVarAnn_idx` (`VarAnnID` ASC),
  CONSTRAINT `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBStudyPar1`
    FOREIGN KEY (`StuParaID`)
    REFERENCES `PAnno`.`StudyParameter` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBVarAnn1`
    FOREIGN KEY (`VarAnnID`)
    REFERENCES `PAnno`.`VarAnn` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB
COMMENT = 'PharmGKBVarAnn to ClinVarStudyParameter';


-- -----------------------------------------------------
-- Table `PAnno`.`VariantSynonyms`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`VariantSynonyms` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Synonym` VARCHAR(256) NOT NULL,
  `VariantID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_VariantSynonyms_Variant1_idx` (`VariantID` ASC),
  CONSTRAINT `fk_VariantSynonyms_Variant1`
    FOREIGN KEY (`VariantID`)
    REFERENCES `PAnno`.`Variant` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`GeneHasVariant`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`GeneHasVariant` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `GeneID` INT(8) UNSIGNED NOT NULL,
  `VariantID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_table1_Gene1_idx` (`GeneID` ASC),
  INDEX `fk_table1_Variant1_idx` (`VariantID` ASC),
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  CONSTRAINT `fk_table1_Gene1`
    FOREIGN KEY (`GeneID`)
    REFERENCES `PAnno`.`Gene` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_table1_Variant1`
    FOREIGN KEY (`VariantID`)
    REFERENCES `PAnno`.`Variant` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`PhenotypeSynonyms`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`PhenotypeSynonyms` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Synonym` VARCHAR(256) NOT NULL,
  `PhenotypeID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_PhenotypeSynonyms_Phenotype1_idx` (`PhenotypeID` ASC),
  CONSTRAINT `fk_PhenotypeSynonyms_Phenotype1`
    FOREIGN KEY (`PhenotypeID`)
    REFERENCES `PAnno`.`Phenotype` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`DrugSynonyms`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`DrugSynonyms` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Synonym` VARCHAR(1024) NOT NULL,
  `DrugID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_DrugSynonyms_Drug1_idx` (`DrugID` ASC),
  CONSTRAINT `fk_DrugSynonyms_Drug1`
    FOREIGN KEY (`DrugID`)
    REFERENCES `PAnno`.`Drug` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`ClinAnnEvidence`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`ClinAnnEvidence` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `CAID` VARCHAR(45) NOT NULL COMMENT 'Clinical Annotation ID',
  `EvidenceID` VARCHAR(64) NOT NULL,
  `TypeID` INT(8) UNSIGNED NOT NULL,
  `URL` VARCHAR(256) NOT NULL,
  `PMID` VARCHAR(32) NULL DEFAULT NULL,
  `Summary` TEXT NULL DEFAULT NULL,
  `Score` TEXT NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnnEvidence_EvidenceTypeDic1_idx` (`TypeID` ASC),
  CONSTRAINT `fk_ClinAnnEvidence_EvidenceTypeDic1`
    FOREIGN KEY (`TypeID`)
    REFERENCES `PAnno`.`EvidenceTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`Guideline`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`Guideline` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Source` VARCHAR(32) NOT NULL,
  `PAID` VARCHAR(16) NOT NULL,
  `Summary` TEXT NOT NULL,
  `Gene` VARCHAR(32) NOT NULL,
  `Drug` VARCHAR(128) NOT NULL,
  `GeneID` INT(8) UNSIGNED NOT NULL,
  `DrugID` INT(8) UNSIGNED NOT NULL,
  `Alternate` TINYINT(1) NOT NULL DEFAULT 0,
  `Dosing` TINYINT(1) NOT NULL DEFAULT 0,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinDosingGuideline_Gene1_idx` (`GeneID` ASC),
  INDEX `fk_ClinDosingGuideline_Drug1_idx` (`DrugID` ASC),
  CONSTRAINT `fk_ClinDosingGuideline_Gene1`
    FOREIGN KEY (`GeneID`)
    REFERENCES `PAnno`.`Gene` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinDosingGuideline_Drug1`
    FOREIGN KEY (`DrugID`)
    REFERENCES `PAnno`.`Drug` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`GuidelineDetail`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`GuidelineDetail` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `PAID` VARCHAR(16) NULL,
  `Phenotype` VARCHAR(256) NOT NULL,
  `Genotype` VARCHAR(256) NULL,
  `Recommendation` TEXT NOT NULL,
  `Avoid` TINYINT(1) NOT NULL DEFAULT 0,
  `Gene` VARCHAR(16) NOT NULL,
  `Drug` VARCHAR(128) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`Domain2Meta`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`Domain2Meta` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `MetaID` INT(8) NOT NULL,
  `DomainID` INT(8) NOT NULL,
  `Class` VARCHAR(32) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`AlleleFunctionality`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`AlleleFunctionality` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Gene` VARCHAR(32) NOT NULL,
  `Allele` VARCHAR(64) NOT NULL,
  `ActivityScore` VARCHAR(32) NULL DEFAULT NULL,
  `Function` VARCHAR(64) NOT NULL,
  `PMID` VARCHAR(512) NULL,
  `Variant` VARCHAR(32) NULL DEFAULT NULL,
  `AlleleManual` VARCHAR(64) NULL DEFAULT NULL,
  `FunctionManual` VARCHAR(64) NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`EHRConsultation`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`EHRConsultation` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'Consultation (Interpretation) Text Provided with Test Result',
  `Gene` VARCHAR(32) NOT NULL,
  `Phenotype` VARCHAR(64) NOT NULL,
  `ActivityScore` VARCHAR(32) NOT NULL,
  `EHR` VARCHAR(64) NOT NULL,
  `Consultation` TEXT NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`DiplotypePhenotype`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`DiplotypePhenotype` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Gene` VARCHAR(32) NOT NULL,
  `Allele1` VARCHAR(64) NOT NULL,
  `Allele2` VARCHAR(64) NOT NULL,
  `ActivityScore` VARCHAR(32) NULL,
  `Phenotype` VARCHAR(64) NOT NULL COMMENT 'Coded Diplotype/Phenotype Summary',
  `EHR` VARCHAR(64) NOT NULL COMMENT 'EHR Priority Notation',
  `ConsultationID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_DiplotypePhenotype_EHRConsultation1_idx` (`ConsultationID` ASC),
  CONSTRAINT `fk_DiplotypePhenotype_EHRConsultation1`
    FOREIGN KEY (`ConsultationID`)
    REFERENCES `PAnno`.`EHRConsultation` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`GuidelineMerge`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`GuidelineMerge` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Source` VARCHAR(32) NOT NULL,
  `PAID` VARCHAR(16) NOT NULL,
  `Summary` TEXT NOT NULL,
  `Phenotype` VARCHAR(256) NOT NULL,
  `Genotype` VARCHAR(256) NULL,
  `Recommendation` TEXT NOT NULL,
  `Avoid` TINYINT(1) NOT NULL DEFAULT 0,
  `Alternate` TINYINT(1) NOT NULL DEFAULT 0,
  `Dosing` TINYINT(1) NOT NULL DEFAULT 0,
  `Gene` VARCHAR(16) NOT NULL,
  `Drug` VARCHAR(128) NOT NULL,
  `GeneID` INT(8) NOT NULL,
  `DrugID` INT(8) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PAnno`.`GuidelineRule`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PAnno`.`GuidelineRule` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Gene` VARCHAR(16) NOT NULL,
  `Variant` VARCHAR(64) NOT NULL,
  `Allele1` VARCHAR(64) NOT NULL,
  `Allele2` VARCHAR(64) NULL,
  `Phenotype` VARCHAR(512) NULL,
  `ClinAnnID` VARCHAR(256) NULL,
  `GuidelineID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_GuidelineRule_GuidelineMerge1_idx` (`GuidelineID` ASC),
  CONSTRAINT `fk_GuidelineRule_GuidelineMerge1`
    FOREIGN KEY (`GuidelineID`)
    REFERENCES `PAnno`.`GuidelineMerge` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
