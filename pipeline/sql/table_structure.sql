-- MySQL dump 10.16  Distrib 10.1.36-MariaDB, for Linux (x86_64)
--
-- Host: localhost    Database: current_project
-- ------------------------------------------------------
-- Server version	10.1.36-MariaDB

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `38clinvar`
--

DROP TABLE IF EXISTS `38clinvar`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `38clinvar` (
  `clinvar_id` int(11) NOT NULL AUTO_INCREMENT,
  `gene_id` int(11) NOT NULL,
  `allele_id` int(11) NOT NULL,
  `symbol` varchar(255) NOT NULL,
  `rs_number` varchar(255) NOT NULL,
  `rcv_accession` varchar(255) NOT NULL,
  `variation_id` int(11) NOT NULL,
  `variant_type` varchar(255) NOT NULL,
  `name` varchar(255) NOT NULL,
  `position` int(11) NOT NULL,
  `variation` varchar(255) NOT NULL,
  `clinical_significance` varchar(255) NOT NULL,
  `last_updated` varchar(255) NOT NULL,
  `phenotypes` varchar(255) NOT NULL,
  `cytogenetic` varchar(255) NOT NULL,
  `review_status` varchar(255) NOT NULL,
  PRIMARY KEY (`clinvar_id`),
  KEY `gene_id` (`gene_id`),
  KEY `gene_id_2` (`gene_id`),
  KEY `gene_id_3` (`gene_id`)
) ENGINE=InnoDB AUTO_INCREMENT=330030 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `CosmicMutantExport`
--

DROP TABLE IF EXISTS `CosmicMutantExport`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `CosmicMutantExport` (
  `cme_id` int(11) NOT NULL AUTO_INCREMENT,
  `gene_name` varchar(255) DEFAULT NULL,
  `accession_number` varchar(255) DEFAULT NULL,
  `gene_cds_length` int(11) DEFAULT NULL,
  `hgnc_id` varchar(255) DEFAULT NULL,
  `sample_name` varchar(255) DEFAULT NULL,
  `id_sample` int(11) DEFAULT NULL,
  `id_tumour` int(11) DEFAULT NULL,
  `primary_site` varchar(255) DEFAULT NULL,
  `site_subtype_1` varchar(255) DEFAULT NULL,
  `site_subtype_2` varchar(255) DEFAULT NULL,
  `site_subtype_3` varchar(255) DEFAULT NULL,
  `primary_histology` varchar(255) DEFAULT NULL,
  `histology_subtype_1` varchar(255) DEFAULT NULL,
  `histology_subtype_2` varchar(255) DEFAULT NULL,
  `histology_subtype_3` varchar(255) DEFAULT NULL,
  `genome_wide_screen` varchar(255) DEFAULT NULL,
  `mutation_id` varchar(255) DEFAULT NULL,
  `mutation_cds` varchar(255) DEFAULT NULL,
  `mutation_aa` varchar(255) DEFAULT NULL,
  `mutation_description` varchar(255) DEFAULT NULL,
  `mutation_zygosity` varchar(255) DEFAULT NULL,
  `loh` varchar(255) DEFAULT NULL,
  `grch` varchar(255) DEFAULT NULL,
  `mutation_genome_position` varchar(255) DEFAULT NULL,
  `mutation_strand` varchar(255) DEFAULT NULL,
  `snp` varchar(255) DEFAULT NULL,
  `resistance_mutation` varchar(255) DEFAULT NULL,
  `fathmm_prediction` varchar(255) DEFAULT NULL,
  `fathmm_score` varchar(255) DEFAULT NULL,
  `mutation_somatic_status` varchar(255) DEFAULT NULL,
  `pubmed_pmid` varchar(255) DEFAULT NULL,
  `id_study` varchar(255) DEFAULT NULL,
  `sample_type` varchar(255) DEFAULT NULL,
  `tumour_origin` varchar(255) DEFAULT NULL,
  `age` varchar(255) DEFAULT NULL,
  `position` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`cme_id`),
  KEY `accession_number` (`accession_number`),
  KEY `gene_name` (`gene_name`),
  KEY `sample_name` (`sample_name`),
  KEY `primary_site` (`primary_site`),
  KEY `primary_histology` (`primary_histology`),
  KEY `mutation_id` (`mutation_id`),
  KEY `mutation_cds` (`mutation_cds`),
  KEY `mutation_aa` (`mutation_aa`),
  KEY `mutation_description` (`mutation_description`),
  KEY `fathmm_prediction` (`fathmm_prediction`),
  KEY `fathmm_score` (`fathmm_score`),
  KEY `mutation_somatic_status` (`mutation_somatic_status`),
  KEY `pubmed_pmid` (`pubmed_pmid`),
  KEY `tumour_origin` (`tumour_origin`),
  KEY `position` (`position`)
) ENGINE=MyISAM AUTO_INCREMENT=6842628 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `clinvar`
--

DROP TABLE IF EXISTS `clinvar`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `clinvar` (
  `clinvar_id` bigint(20) DEFAULT NULL,
  `gene_id` bigint(20) DEFAULT NULL,
  `allele_id` bigint(20) DEFAULT NULL,
  `symbol` varchar(255) DEFAULT NULL,
  `rs_number` varchar(255) DEFAULT NULL,
  `rcv_accession` varchar(255) DEFAULT NULL,
  `variation_id` bigint(20) DEFAULT NULL,
  `variant_type` varchar(255) DEFAULT NULL,
  `name` varchar(255) DEFAULT NULL,
  `clinical_significance` varchar(255) DEFAULT NULL,
  `last_updated` varchar(255) DEFAULT NULL,
  `cytogenetic` varchar(255) DEFAULT NULL,
  `review_status` varchar(255) DEFAULT NULL,
  `phenotypes` varchar(255) DEFAULT NULL,
  `nm_id` text,
  `variation` varchar(255) DEFAULT NULL,
  `position` varchar(255) DEFAULT '0',
  `np_id` varchar(255) DEFAULT NULL,
  KEY `ix_clinvar_clinvar_id` (`clinvar_id`),
  KEY `allele_id` (`allele_id`),
  KEY `allele_id_2` (`allele_id`),
  KEY `symbol` (`symbol`),
  KEY `rs_number` (`rs_number`),
  KEY `rcv_accession` (`rcv_accession`),
  KEY `variation_id` (`variation_id`),
  KEY `variant_type` (`variant_type`),
  KEY `variant_type_2` (`variant_type`),
  KEY `name` (`name`),
  KEY `name_2` (`name`),
  KEY `clinical_significance` (`clinical_significance`),
  KEY `clinical_significance_2` (`clinical_significance`),
  KEY `last_updated` (`last_updated`),
  KEY `last_updated_2` (`last_updated`),
  KEY `last_updated_3` (`last_updated`),
  KEY `cytogenetic` (`cytogenetic`),
  KEY `cytogenetic_2` (`cytogenetic`),
  KEY `cytogenetic_3` (`cytogenetic`),
  KEY `review_status` (`review_status`),
  KEY `review_status_2` (`review_status`),
  KEY `review_status_3` (`review_status`),
  KEY `review_status_4` (`review_status`),
  KEY `phenotypes` (`phenotypes`),
  KEY `phenotypes_2` (`phenotypes`),
  KEY `phenotypes_3` (`phenotypes`),
  KEY `phenotypes_4` (`phenotypes`),
  KEY `variation` (`variation`),
  KEY `variation_2` (`variation`),
  KEY `variation_3` (`variation`),
  KEY `variation_4` (`variation`),
  KEY `position` (`position`),
  KEY `position_2` (`position`),
  KEY `position_3` (`position`),
  KEY `position_4` (`position`),
  KEY `np_id` (`np_id`),
  KEY `np_id_2` (`np_id`),
  KEY `np_id_3` (`np_id`),
  KEY `np_id_4` (`np_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `conservation_scores`
--

DROP TABLE IF EXISTS `conservation_scores`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `conservation_scores` (
  `cs_id` int(11) NOT NULL AUTO_INCREMENT,
  `msa_id` int(11) NOT NULL,
  `aminoacid_number` int(11) NOT NULL,
  `specie` varchar(255) NOT NULL,
  `transcript_id_specie` varchar(255) NOT NULL,
  `score` int(11) NOT NULL,
  `score_type` varchar(255) NOT NULL,
  `transcript_id` varchar(255) NOT NULL,
  `gene_symbol` varchar(255) NOT NULL,
  PRIMARY KEY (`cs_id`),
  KEY `transcript_id` (`transcript_id`),
  KEY `msa_id` (`msa_id`)
) ENGINE=MyISAM AUTO_INCREMENT=2062333 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `convart_gene`
--

DROP TABLE IF EXISTS `convart_gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `convart_gene` (
  `id` int(32) NOT NULL,
  `sequence` text NOT NULL,
  `species_id` varchar(55) NOT NULL,
  `hash` varchar(150) NOT NULL
) ENGINE=MyISAM DEFAULT CHARSET=utf8;

--
-- Indexes for table `convart_gene`
--
ALTER TABLE `convart_gene`
  ADD PRIMARY KEY (`id`),
  ADD UNIQUE KEY `hash` (`hash`),
  ADD KEY `species` (`species_id`);

/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `convart_gene_to_db`
--

DROP TABLE IF EXISTS `convart_gene_to_db`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `convart_gene_to_db` (
  `convart_gene_id` int(11) NOT NULL,
  `db` varchar(25) NOT NULL,
  `db_id` varchar(50) NOT NULL,
  `db_id_version` int(11) DEFAULT NULL,
  PRIMARY KEY (`convart_gene_id`,`db`,`db_id`),
  KEY `db` (`db`),
  KEY `db_id` (`db_id`),
  KEY `convart_gene_id` (`convart_gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `dbsnp`
--

DROP TABLE IF EXISTS `dbsnp`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `dbsnp` (
  `dbsnp_index` bigint(20) DEFAULT NULL,
  `Uploaded_variation` text,
  `Location` text,
  `Allele` text,
  `Gene` text,
  `Feature` varchar(255) DEFAULT NULL,
  `Feature_type` text,
  `Consequence` text,
  `cDNA_position` text,
  `CDS_position` text,
  `Protein_position` text,
  `HGVSp` text,
  `Impact` text,
  KEY `ix_dbsnp_index` (`dbsnp_index`),
  KEY `Feature` (`Feature`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `disease_info`
--

DROP TABLE IF EXISTS `disease_info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `disease_info` (
  `di_id` int(11) NOT NULL AUTO_INCREMENT,
  `disease_id` varchar(255) NOT NULL,
  `disease_name` varchar(255) NOT NULL,
  `category` varchar(255) NOT NULL,
  PRIMARY KEY (`di_id`),
  UNIQUE KEY `disease_id` (`disease_id`),
  KEY `disease_name` (`disease_name`)
) ENGINE=InnoDB AUTO_INCREMENT=10371 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `diseases_genes`
--

DROP TABLE IF EXISTS `diseases_genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `diseases_genes` (
  `dg_id` int(11) NOT NULL AUTO_INCREMENT,
  `gene_id` int(11) NOT NULL,
  `disease_id` varchar(255) NOT NULL,
  `dsi` float NOT NULL,
  `dpi` float NOT NULL,
  PRIMARY KEY (`dg_id`),
  KEY `gene_id` (`gene_id`)
) ENGINE=InnoDB AUTO_INCREMENT=81747 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `domains`
--

DROP TABLE IF EXISTS `domains`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `domains` (
  `domain_id` int(11) NOT NULL AUTO_INCREMENT,
  `transcript_id` varchar(255) NOT NULL,
  `pfam_id` varchar(255) NOT NULL,
  `pfam_name` varchar(255) NOT NULL,
  `type` varchar(255) NOT NULL,
  `start_point` int(11) NOT NULL,
  `end_point` int(11) NOT NULL,
  `bitscore` float NOT NULL,
  `evalue` varchar(10) NOT NULL,
  `clan` varchar(255) NOT NULL,
  PRIMARY KEY (`domain_id`),
  KEY `transcript_id` (`transcript_id`)
) ENGINE=InnoDB AUTO_INCREMENT=131071 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `domains_desc`
--

DROP TABLE IF EXISTS `domains_desc`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `domains_desc` (
  `domain_desc_id` int(11) NOT NULL AUTO_INCREMENT,
  `pfam_id` varchar(11) NOT NULL,
  `pfam_name` varchar(255) NOT NULL,
  `domain_desc` varchar(255) NOT NULL,
  PRIMARY KEY (`domain_desc_id`),
  UNIQUE KEY `pfam_id` (`pfam_id`)
) ENGINE=InnoDB AUTO_INCREMENT=17930 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gnomad`
--

DROP TABLE IF EXISTS `gnomad`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gnomad` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `consequence` varchar(255) NOT NULL,
  `ensg_id` varchar(255) NOT NULL,
  `canonical_transcript` varchar(80) NOT NULL,
  `variation` varchar(50) NOT NULL,
  `filters` varchar(255) NOT NULL,
  `rs_number` varchar(50) NOT NULL,
  `variant_id` varchar(255) CHARACTER SET latin1 COLLATE latin1_general_ci NOT NULL,
  `position` int(11) NOT NULL,
  `allele_count` int(11) NOT NULL,
  `allele_number` int(11) NOT NULL,
  `allelle_frequency` int(11) NOT NULL,
  `flags` varchar(255) NOT NULL,
  `datasets` varchar(512) NOT NULL,
  `is_canon` tinyint(4) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `canonical_transcript_2` (`canonical_transcript`,`variation`,`rs_number`),
  KEY `canonical_transcript` (`canonical_transcript`),
  KEY `ensg_id` (`ensg_id`),
  KEY `rs_number` (`rs_number`),
  KEY `consequence` (`consequence`),
  KEY `is_canon` (`is_canon`)
) ENGINE=MyISAM AUTO_INCREMENT=17247788 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `homology`
--

DROP TABLE IF EXISTS `homology`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `homology` (
  `homology_id` int(11) NOT NULL AUTO_INCREMENT,
  `human_gene_id` int(11) NOT NULL,
  `chimp_gene_id` varchar(255) NOT NULL,
  `macaque_gene_id` varchar(255) NOT NULL,
  `rat_gene_id` varchar(255) NOT NULL,
  `mouse_gene_id` varchar(255) NOT NULL,
  `zebrafish_gene_id` varchar(255) NOT NULL,
  `frog_gene_id` varchar(255) NOT NULL,
  `fruitfly_gene_id` varchar(255) NOT NULL,
  `worm_gene_id` varchar(255) NOT NULL,
  `yeast_gene_id` varchar(255) NOT NULL,
  PRIMARY KEY (`homology_id`)
) ENGINE=MyISAM AUTO_INCREMENT=32768 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mapping_human`
--

DROP TABLE IF EXISTS `mapping_human`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mapping_human` (
  `mh_id` int(11) NOT NULL AUTO_INCREMENT,
  `gene_id` int(11) NOT NULL,
  `other_ids` varchar(255) NOT NULL,
  `gene_symbol` varchar(255) NOT NULL,
  `gene_synonyms` varchar(255) NOT NULL,
  `protein_numbers` varchar(255) NOT NULL,
  `gene_description` varchar(255) NOT NULL,
  PRIMARY KEY (`mh_id`),
  UNIQUE KEY `gene_id` (`gene_id`)
) ENGINE=InnoDB AUTO_INCREMENT=19264 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `msa`
--

DROP TABLE IF EXISTS `msa`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `msa` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `fasta` mediumtext COLLATE latin1_bin NOT NULL,
  `alignment_method` varchar(255) CHARACTER SET utf8 NOT NULL,
  `created_date` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=MyISAM AUTO_INCREMENT=55543 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `msa_best_combination`
--

DROP TABLE IF EXISTS `msa_best_combination`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `msa_best_combination` (
  `msa_id` int(11) NOT NULL,
  `convart_gene_id` int(11) NOT NULL,
  PRIMARY KEY (`msa_id`),
  KEY `convart_gene_id` (`convart_gene_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `msa_gene`
--

DROP TABLE IF EXISTS `msa_gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `msa_gene` (
  `msa_id` int(11) NOT NULL,
  `convart_gene_id` int(11) NOT NULL,
  PRIMARY KEY (`msa_id`,`convart_gene_id`),
  KEY `convart_gene_id` (`convart_gene_id`),
  KEY `msa_id` (`msa_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ncbi_gene_meta`
--

DROP TABLE IF EXISTS `ncbi_gene_meta`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ncbi_gene_meta` (
  `ncbi_gene_id` int(11) NOT NULL,
  `meta_key` varchar(50) NOT NULL,
  `meta_value` varchar(255) NOT NULL,
  PRIMARY KEY (`ncbi_gene_id`,`meta_key`,`meta_value`),
  KEY `meta_key` (`meta_key`),
  KEY `meta_value` (`meta_value`) USING HASH,
  KEY `ncbi_gene_id_2` (`ncbi_gene_id`) USING HASH
) ENGINE=MyISAM DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ptm`
--

DROP TABLE IF EXISTS `ptm`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ptm` (
  `index_id` bigint(20) DEFAULT NULL,
  `gene` text,
  `protein` text,
  `acc_id` varchar(32) DEFAULT NULL,
  `hu_chr_loc` text,
  `mod_rsd` text,
  `site_grp_id` bigint(20) DEFAULT NULL,
  `organism` text,
  `mw_kd` double DEFAULT NULL,
  `domain` text,
  `site_+/-7_aa` text,
  `lt_lit` double DEFAULT NULL,
  `ms_lit` double DEFAULT NULL,
  `ms_cst` double DEFAULT NULL,
  `cst_cat#` text,
  `ptm_type` text,
  `enst_id` text,
  `position` text,
  KEY `ix_ptm_index` (`index_id`),
  KEY `acc_id` (`acc_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `species` (
  `id` varchar(255) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2019-08-05 14:10:42
