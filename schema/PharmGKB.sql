-- MySQL Script generated by MySQL Workbench
-- Thu Dec  9 19:20:34 2021
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';


-- -----------------------------------------------------
-- Schema PharmGKB
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `PharmGKB` DEFAULT CHARACTER SET utf8 ;
USE `PharmGKB` ;

-- -----------------------------------------------------
-- Table `PharmGKB`.`CharacteristicsTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`CharacteristicsTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(32) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB
COMMENT = 'Charicteristic type dictionary';


-- -----------------------------------------------------
-- Table `PharmGKB`.`Gene`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`Gene` (
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
-- Table `PharmGKB`.`Variant`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`Variant` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL COMMENT 'Variant(rsID) or haplotype',
  `NCID` VARCHAR(32) NULL DEFAULT NULL,
  `PAID` VARCHAR(16) NULL DEFAULT NULL,
  `Type` VARCHAR(16) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`Drug`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`Drug` (
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
-- Table `PharmGKB`.`PhenotypeCategoryDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`PhenotypeCategoryDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`Phenotype`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`Phenotype` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(128) NOT NULL,
  `PAID` VARCHAR(16) NOT NULL,
  `Type` VARCHAR(32) NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `PAID_UNIQUE` (`PAID` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`ClinAnn`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`ClinAnn` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT COMMENT 'ClinicalAnn',
  `CAID` INT(16) UNSIGNED NOT NULL COMMENT 'Clinical Annotation ID',
  `Gene` VARCHAR(64) NOT NULL,
  `VariantOrHaplotype` VARCHAR(256) NOT NULL,
  `Drug` VARCHAR(256) NOT NULL,
  `Phenotypes` VARCHAR(256) NOT NULL,
  `EvidenceLevel` VARCHAR(4) NULL DEFAULT NULL,
  `LevelOverride` TEXT NULL DEFAULT NULL,
  `LevelModifier` VARCHAR(32) NULL DEFAULT NULL,
  `Score` VARCHAR(32) NULL DEFAULT NULL,
  `PhenotypeCategoryID` INT(8) UNSIGNED NOT NULL,
  `GenotypeOrAllele` VARCHAR(256) NOT NULL,
  `Annotation` TEXT NOT NULL,
  `Function` VARCHAR(64) NOT NULL,
  `PMIDCount` INT(4) NOT NULL,
  `EvidenceCount` INT(4) NOT NULL,
  `URL` VARCHAR(256) NOT NULL,
  `SpecialtyPopulation` VARCHAR(64) NULL DEFAULT NULL,
  `GeneID` INT(8) UNSIGNED NOT NULL,
  `VariantID` INT(8) UNSIGNED NOT NULL,
  `DrugID` INT(8) UNSIGNED NOT NULL,
  `Comment` TEXT NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnn_Gene1_idx` (`GeneID` ASC),
  INDEX `fk_ClinAnn_Variant1_idx` (`VariantID` ASC),
  INDEX `fk_ClinAnn_Drug1_idx` (`DrugID` ASC),
  INDEX `fk_ClinAnn_PhenotypeCategoryDic1_idx` (`PhenotypeCategoryID` ASC),
  CONSTRAINT `fk_ClinAnn_Gene1`
    FOREIGN KEY (`GeneID`)
    REFERENCES `PharmGKB`.`Gene` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnn_Variant1`
    FOREIGN KEY (`VariantID`)
    REFERENCES `PharmGKB`.`Variant` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnn_Drug1`
    FOREIGN KEY (`DrugID`)
    REFERENCES `PharmGKB`.`Drug` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnn_PhenotypeCategoryDic1`
    FOREIGN KEY (`PhenotypeCategoryID`)
    REFERENCES `PharmGKB`.`PhenotypeCategoryDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`RaceDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`RaceDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB
COMMENT = 'dictionary: races';


-- -----------------------------------------------------
-- Table `PharmGKB`.`EvidenceTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`EvidenceTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`VarAnn`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`VarAnn` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `AnnotationID` INT(16) UNSIGNED NOT NULL,
  `VariantOrHaplotype` VARCHAR(512) NOT NULL,
  `Gene` VARCHAR(512) NOT NULL,
  `Drug` VARCHAR(512) NOT NULL,
  `PMID` VARCHAR(64) NULL DEFAULT NULL COMMENT 'PMID of article describing the relation',
  `PhenotypeCategoryID` INT(8) UNSIGNED NOT NULL,
  `Significance` VARCHAR(16) NOT NULL COMMENT 'Foreign key -> PharmGKBSignificance.ID',
  `Note` TEXT NULL DEFAULT NULL COMMENT 'Full describtion of the relation in different heterozygotes/homozygotes',
  `Sentence` TEXT NULL DEFAULT NULL COMMENT 'Summary of Notes',
  `Allele` VARCHAR(256) NULL DEFAULT NULL COMMENT 'Allele studied in the article',
  `SpecialtyPopulation` VARCHAR(64) NULL DEFAULT NULL,
  `TypeID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  FULLTEXT INDEX `Notes` (`Note`),
  FULLTEXT INDEX `Sentence` (`Sentence`),
  INDEX `fk_VarAnn_EvidenceTypeDic1_idx` (`TypeID` ASC),
  INDEX `fk_VarAnn_PhenotypeCategoryDic1_idx` (`PhenotypeCategoryID` ASC),
  CONSTRAINT `fk_VarAnn_EvidenceTypeDic1`
    FOREIGN KEY (`TypeID`)
    REFERENCES `PharmGKB`.`EvidenceTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_VarAnn_PhenotypeCategoryDic1`
    FOREIGN KEY (`PhenotypeCategoryID`)
    REFERENCES `PharmGKB`.`PhenotypeCategoryDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`RatioStatTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`RatioStatTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(8) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`StudyParameter`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`StudyParameter` (
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
    REFERENCES `PharmGKB`.`CharacteristicsTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBStudyParameter_PharmGKBRatioStatType1`
    FOREIGN KEY (`RatioStatTypeID`)
    REFERENCES `PharmGKB`.`RatioStatTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`StuParaHasRace`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`StuParaHasRace` (
  `StuParaID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `RaceID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBRace1_idx` (`RaceID` ASC),
  INDEX `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBStudyPar_idx` (`StuParaID` ASC),
  CONSTRAINT `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBRace1`
    FOREIGN KEY (`RaceID`)
    REFERENCES `PharmGKB`.`RaceDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBStudyParameter_has_PharmGKBRace_PharmGKBStudyParam1`
    FOREIGN KEY (`StuParaID`)
    REFERENCES `PharmGKB`.`StudyParameter` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`StudyTypeDic`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`StudyTypeDic` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(64) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  UNIQUE INDEX `Name_UNIQUE` (`Name` ASC))
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`StuParaHasStuType`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`StuParaHasStuType` (
  `StuTypeID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `StuParaID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStu_idx` (`StuParaID` ASC),
  INDEX `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStu_idx1` (`StuTypeID` ASC),
  CONSTRAINT `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStudy1`
    FOREIGN KEY (`StuTypeID`)
    REFERENCES `PharmGKB`.`StudyTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBStudyType_has_PharmGKBStudyParameter_PharmGKBStudy2`
    FOREIGN KEY (`StuParaID`)
    REFERENCES `PharmGKB`.`StudyParameter` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`VarAnnHasStuPara`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`VarAnnHasStuPara` (
  `VarAnnID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `StuParaID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBStudyP_idx` (`StuParaID` ASC),
  INDEX `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBVarAnn_idx` (`VarAnnID` ASC),
  CONSTRAINT `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBStudyPar1`
    FOREIGN KEY (`StuParaID`)
    REFERENCES `PharmGKB`.`StudyParameter` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_PharmGKBVarAnn_has_PharmGKBStudyParameter_PharmGKBVarAnn1`
    FOREIGN KEY (`VarAnnID`)
    REFERENCES `PharmGKB`.`VarAnn` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`VariantSynonyms`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`VariantSynonyms` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Synonym` VARCHAR(256) NOT NULL,
  `VariantID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_VariantSynonyms_Variant1_idx` (`VariantID` ASC),
  CONSTRAINT `fk_VariantSynonyms_Variant1`
    FOREIGN KEY (`VariantID`)
    REFERENCES `PharmGKB`.`Variant` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`GeneHasVariant`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`GeneHasVariant` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `GeneID` INT(8) UNSIGNED NOT NULL,
  `VariantID` INT(8) UNSIGNED NOT NULL,
  INDEX `fk_table1_Gene1_idx` (`GeneID` ASC),
  INDEX `fk_table1_Variant1_idx` (`VariantID` ASC),
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  CONSTRAINT `fk_table1_Gene1`
    FOREIGN KEY (`GeneID`)
    REFERENCES `PharmGKB`.`Gene` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_table1_Variant1`
    FOREIGN KEY (`VariantID`)
    REFERENCES `PharmGKB`.`Variant` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`PhenotypeSynonyms`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`PhenotypeSynonyms` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Synonym` VARCHAR(256) NOT NULL,
  `PhenotypeID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_PhenotypeSynonyms_Phenotype1_idx` (`PhenotypeID` ASC),
  CONSTRAINT `fk_PhenotypeSynonyms_Phenotype1`
    FOREIGN KEY (`PhenotypeID`)
    REFERENCES `PharmGKB`.`Phenotype` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`DrugSynonyms`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`DrugSynonyms` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Synonym` VARCHAR(1024) NOT NULL,
  `DrugID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_DrugSynonyms_Drug1_idx` (`DrugID` ASC),
  CONSTRAINT `fk_DrugSynonyms_Drug1`
    FOREIGN KEY (`DrugID`)
    REFERENCES `PharmGKB`.`Drug` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`ClinAnnEvidence`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`ClinAnnEvidence` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `CAID` VARCHAR(45) NOT NULL COMMENT 'Clinical Annotation ID',
  `EvidenceID` VARCHAR(64) NOT NULL,
  `EvidenceTypeID` INT(8) UNSIGNED NOT NULL,
  `URL` VARCHAR(256) NOT NULL,
  `PMID` VARCHAR(32) NULL DEFAULT NULL,
  `Summary` TEXT NULL DEFAULT NULL,
  `Score` TEXT NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnnEvidence_EvidenceTypeDic1_idx` (`EvidenceTypeID` ASC),
  CONSTRAINT `fk_ClinAnnEvidence_EvidenceTypeDic1`
    FOREIGN KEY (`EvidenceTypeID`)
    REFERENCES `PharmGKB`.`EvidenceTypeDic` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`ClinAnnHasPhenotype`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`ClinAnnHasPhenotype` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `ClinAnnID` INT(8) UNSIGNED NOT NULL,
  `PhenotypeID` INT(8) UNSIGNED NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinAnnHasPhenotype_ClinAnn1_idx` (`ClinAnnID` ASC),
  INDEX `fk_ClinAnnHasPhenotype_Phenotype1_idx` (`PhenotypeID` ASC),
  CONSTRAINT `fk_ClinAnnHasPhenotype_ClinAnn1`
    FOREIGN KEY (`ClinAnnID`)
    REFERENCES `PharmGKB`.`ClinAnn` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinAnnHasPhenotype_Phenotype1`
    FOREIGN KEY (`PhenotypeID`)
    REFERENCES `PharmGKB`.`Phenotype` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`ClinDosingGuideline`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`ClinDosingGuideline` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `Source` VARCHAR(32) NOT NULL,
  `Annotation` TEXT NOT NULL,
  `RelatedGeneID` INT(8) UNSIGNED NOT NULL,
  `RelatedDrugID` INT(8) UNSIGNED NOT NULL,
  `URL` VARCHAR(256) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC),
  INDEX `fk_ClinDosingGuideline_Gene1_idx` (`RelatedGeneID` ASC),
  INDEX `fk_ClinDosingGuideline_Drug1_idx` (`RelatedDrugID` ASC),
  CONSTRAINT `fk_ClinDosingGuideline_Gene1`
    FOREIGN KEY (`RelatedGeneID`)
    REFERENCES `PharmGKB`.`Gene` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE,
  CONSTRAINT `fk_ClinDosingGuideline_Drug1`
    FOREIGN KEY (`RelatedDrugID`)
    REFERENCES `PharmGKB`.`Drug` (`ID`)
    ON DELETE CASCADE
    ON UPDATE CASCADE)
ENGINE = InnoDB;


-- -----------------------------------------------------
-- Table `PharmGKB`.`Domain2Meta`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `PharmGKB`.`Domain2Meta` (
  `ID` INT(8) UNSIGNED NOT NULL AUTO_INCREMENT,
  `MetaID` INT(8) NOT NULL,
  `DomainID` INT(8) NOT NULL,
  `Class` VARCHAR(32) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC))
ENGINE = InnoDB;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;
