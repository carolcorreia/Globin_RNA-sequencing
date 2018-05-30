######################################
# Globin-derived transcripts RNA-seq #
#       2018 Analysis - Part A       #
######################################

# Author: Carolina N. Correia
# GitHub Repository DOI:
# Date: May 29th 2018

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Define variables for specific directories
salmonDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon"
horseDir <- file.path(paste0(salmonDir, "/TPM/horse_2018_TPM"))
cattleDir <- file.path(paste0(salmonDir, "/TPM/cattle_2018_TPM"))

gffDir <- file.path(paste0(getwd(), "/GFF3_RefSeq"))
horseGFF <- file.path(paste0(gffDir, "/ref_EquCab3.0_top_level.gff3"))
cattleGFF <- file.path(paste0(gffDir, "/ref_ARS-UCD1.2_top_level.gff3"))

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(plyr)
library(tidyverse)
library(devtools)
library(magrittr)
library(stringr)
library(forcats)
library(rtracklayer)
library(AnnotationDbi)
library(GenomicFeatures)
library(tximport)
library(rjson)
library(reshape2)


# Uncomment functions below to install packages in case you don't have them

# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#biocLite("tximport")
#biocLite("GenomicFeatures")


# CRAN packages
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("rjson")
#install.packages("reshape2")

################################################################
# 03 Create NCBI RefSeq TxDb objects for transcript annotation #
#                      (horse and cow)                         #
################################################################

# Cow NCBI RefSeq .gff3 file used in this analysis
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Bos_taurus/GFF/ref_ARS-UCD1.2_top_level.gff3.gz
# ANNOTATION RELEASE NAME: NCBI Bos taurus Annotation Release 106
# ANNOTATION EVIDENCE FREEZE DATE: 19 April 2018
# ANNOTATION RELEASE DATE: 11 May 2018
# ANNOTATION REPORT: https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Bos_taurus/106/

# Horse NCBI RefSeq .gff3 file used in this analysis
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Equus_caballus/GFF/ref_EquCab3.0_top_level.gff3.gz
# ANNOTATION RELEASE NAME: NCBI Equus caballus Annotation Release 103
# ANNOTATION EVIDENCE FREEZE DATE: 11 January 2018
# ANNOTATION RELEASE DATE: 26 January 2018
# ANNOTATION REPORT: https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Equus_caballus/103/

# Create DBs using .gff3 files
horseDB <- makeTxDbFromGFF(horseGFF, format = "gff3")
cattleDB <- makeTxDbFromGFF(cattleGFF, format = "gff3")

# Check databases info
horseDB
cattleDB

################################################################
# 04 Create dataframe of transcripts and gene IDs for tximport #
################################################################

# Check DB key type
keytypes(horseDB)
keytypes(cattleDB)

# Extract keys from DBs
horse_keys <- keys(horseDB, keytype = "GENEID")
cattle_keys <- keys(cattleDB, keytype = "GENEID")

# Create transcript-to-gene dataframe for each species
horse_tx2gene <- AnnotationDbi::select(horseDB,
                                       keys = horse_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

cattle_tx2gene <- AnnotationDbi::select(cattleDB,
                                        keys = cattle_keys,
                                        keytype = "GENEID",
                                        columns = "TXNAME")

# tximport requires the order to be: transcript name followed by gene ID
horse_tx2gene %<>%
    dplyr::select(TXNAME, GENEID)
head(horse_tx2gene)


cattle_tx2gene %<>%
    dplyr::select(TXNAME, GENEID)
head(cattle_tx2gene)

##############################################################
# 05 Import salmon TPM estimates and summarise at gene-level #
##############################################################

# Get paths to salmon files
horse_files <- list.files(horseDir, pattern = "quant.sf", full.names = TRUE)
names(horse_files) <- list.files(horseDir)

cattle_files <- list.files(cattleDir, pattern = "quant.sf", full.names = TRUE)
names(cattle_files) <- list.files(cattleDir)

# Import transcript-level estimates summarized to the gene-level
# (it's the default when txOut = FALSE or not provided)
horse_txi <- tximport(horse_files,
                      type = "salmon",
                      tx2gene = horse_tx2gene)
names(horse_txi)
head(horse_txi$abundance)


cattle_txi <- tximport(cattle_files,
                       type = "salmon",
                       tx2gene = cattle_tx2gene)
names(cattle_txi)
head(cattle_txi$abundance)

############################################
# 06 Remove zero and lowly expressed genes #
############################################

# Convert gene-level TPM abundances into data frames
horse_TPM <- as.data.frame(horse_txi$abundance)
cattle_TPM <- as.data.frame(cattle_txi$abundance)

# Remove non-expressed genes
horse_nozeros <- horse_TPM[rowSums(horse_TPM) > 0, ]
dim(horse_nozeros)
dim(horse_TPM)

cattle_nozeros <- cattle_TPM[rowSums(cattle_TPM) > 0, ]
dim(cattle_nozeros)
dim(cattle_TPM)

# Remove lowly expressed genes (< 1 TPM in approximately half of total samples)
horse_filt <- horse_nozeros[rowSums(horse_nozeros >= 1) >= 18, ]
dim(horse_filt)
dim(horse_nozeros)

cattle_filt <- cattle_nozeros[rowSums(cattle_nozeros >= 1) >= 5, ]
dim(cattle_filt)
dim(cattle_nozeros)

####################################
# 07 Tidy filtered TPM data frames #
####################################
horse_filt %<>%
    rownames_to_column(var = "Gene_RefSeqID") %>%
    melt(value.name = "TPM") %>%
    dplyr::rename(sample = variable) %>%
    dplyr::mutate(species = "Equine")

cattle_filt %<>%
    rownames_to_column(var = "Gene_RefSeqID") %>%
    melt(value.name = "TPM") %>%
    dplyr::rename(sample = variable) %>%
    dplyr::mutate(species = "Bovine")

# Check reformatted data frames
head(horse_filt)
head(cattle_filt)

#########################################################
# 08 Row-bind filtered TPM data frames from all species #
#########################################################

horse_filt %>%
dplyr::bind_rows(cattle_filt) %>%
    as.tibble() -> TPM_filt_all

# Check total number of rows
length(TPM_filt_all$Gene_RefSeqID) ==
    length(horse_filt$Gene_RefSeqID) +
    length(cattle_filt$Gene_RefSeqID)

# Visualise TPM df
TPM_filt_all

#######################################
# 09 Convert species column to factor #
#######################################

# Order of the levels match the order of the first appearance in the data
TPM_filt_all$species %<>%
    factor() %>%
    forcats::fct_inorder()

# Check species factors
levels(TPM_filt_all$species)

#####################################################
# 10 Add plotting labels to filtered TPM data frame #
#####################################################

# Create new column
TPM_filt_all$labels <- TPM_filt_all$sample

# Correct plotting labels
TPM_filt_all$labels %<>%
    str_replace("_quant.sf", "") %>%
    str_replace("A", "") %>%
    str_replace("_W-1_F", "") %>%
    str_replace("6511", "Bta_01_U") %>%
    str_replace("6514", "Bta_02_U") %>%
    str_replace("6520", "Bta_03_U") %>%
    str_replace("6522", "Bta_04_U") %>%
    str_replace("6526", "Bta_05_U") %>%
    str_replace("6635", "Bta_06_U") %>%
    str_replace("6636", "Bta_07_U") %>%
    str_replace("6637", "Bta_08_U") %>%
    str_replace("6644", "Bta_09_U") %>%
    str_replace("6698", "Bta_10_U") %>%
    str_replace("SRR3671009", "Eca_T01_U") %>%
    str_replace("SRR3671010", "Eca_T02_U") %>%
    str_replace("SRR3671011", "Eca_T03_U") %>%
    str_replace("SRR3671012", "Eca_T04_U") %>%
    str_replace("SRR3671013", "Eca_T05_U") %>%
    str_replace("SRR3671014", "Eca_T06_U") %>%
    str_replace("SRR3671015", "Eca_T07_U") %>%
    str_replace("SRR3671016", "Eca_T08_U") %>%
    str_replace("SRR3671017", "Eca_T09_U") %>%
    str_replace("SRR3671018", "Eca_T10_U") %>%
    str_replace("SRR3671019", "Eca_T11_U") %>%
    str_replace("SRR3671020", "Eca_T12_U") %>%
    str_replace("SRR3671021", "Eca_S13_U") %>%
    str_replace("SRR3671022", "Eca_T14_U") %>%
    str_replace("SRR3671023", "Eca_T15_U") %>%
    str_replace("SRR3671024", "Eca_T16_U") %>%
    str_replace("SRR3671025", "Eca_T17_U") %>%
    str_replace("SRR3671026", "Eca_T18_U") %>%
    str_replace("SRR3671027", "Eca_T19_U") %>%
    str_replace("SRR3671028", "Eca_T20_U") %>%
    str_replace("SRR3671029", "Eca_T21_U") %>%
    str_replace("SRR3671030", "Eca_T22_U") %>%
    str_replace("SRR3671031", "Eca_T23_U") %>%
    str_replace("SRR3671032", "Eca_T24_U") %>%
    str_replace("SRR3671033", "Eca_T25_U") %>%
    str_replace("SRR3671034", "Eca_T26_U") %>%
    str_replace("SRR3671035", "Eca_T27_U") %>%
    str_replace("SRR3671036", "Eca_T28_U") %>%
    str_replace("SRR3671037", "Eca_T29_U") %>%
    str_replace("SRR3671038", "Eca_T30_U") %>%
    str_replace("SRR3671039", "Eca_T31_U") %>%
    str_replace("SRR3671040", "Eca_T32_U") %>%
    str_replace("SRR3671041", "Eca_S33_U") %>%
    str_replace("SRR3671042", "Eca_S34_U") %>%
    str_replace("SRR3671043", "Eca_S35_U") %>%
    str_replace("SRR3671044", "Eca_S36_U") %>%
    str_replace("SRR3671045", "Eca_S37_U")

# Check labels
unique(TPM_filt_all$labels)

# Convert labels into factors and order them
TPM_filt_all$labels %<>%
    factor(levels = c("Eca_T01_U", "Eca_T02_U", "Eca_T03_U", "Eca_T04_U",
                      "Eca_T05_U", "Eca_T06_U", "Eca_T07_U", "Eca_T08_U",
                      "Eca_T09_U", "Eca_T10_U", "Eca_T11_U", "Eca_T12_U",
                      "Eca_S13_U", "Eca_T14_U", "Eca_T15_U", "Eca_T16_U",
                      "Eca_T17_U", "Eca_T18_U", "Eca_T19_U", "Eca_T20_U",
                      "Eca_T21_U", "Eca_T22_U", "Eca_T23_U", "Eca_T24_U",
                      "Eca_T25_U", "Eca_T26_U", "Eca_T27_U", "Eca_T28_U",
                      "Eca_T29_U", "Eca_T30_U", "Eca_T31_U", "Eca_T32_U",
                      "Eca_S33_U", "Eca_S34_U", "Eca_S35_U", "Eca_S36_U",
                      "Eca_S37_U", "Bta_01_U", "Bta_02_U", "Bta_03_U",
                      "Bta_04_U", "Bta_05_U", "Bta_06_U", "Bta_07_U",
                      "Bta_08_U", "Bta_09_U", "Bta_10_U"))

# Check factor order
levels(TPM_filt_all$labels)

####################################################
# 11 Add treatment info to filtered TPM data frame #
####################################################

# Create new column
TPM_filt_all$treatment <- TPM_filt_all$labels

# Correct info and convert into factor
TPM_filt_all$treatment %<>%
    stringr::str_replace("Bta_\\d\\d_", "") %>%
    stringr::str_replace("Eca_(T|S)\\d\\d_", "") %>%
    stringr::str_replace("U", "Undepleted") %>% 
    factor()

# Check treatment factors
levels(TPM_filt_all$treatment)

#######################
# 12 Save .RData file #
#######################

save.image(file = "Globin-2018_analysis.RData")

#########################
# 13 Get R session info #
#########################

devtools::session_info()

######################################
# Proceed to Part B of 2018 analysis #
######################################

# File: 10-Globin-2018_analysis-Part_B.R




















