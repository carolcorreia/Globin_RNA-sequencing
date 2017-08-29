######################################
# Globin-derived transcripts RNA-seq #
#   Additional Analysis - Part A     #
######################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: August 29th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Define variables for specific directories
TPMdir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/extra_horse_TPM"
gffDir <- file.path(paste0(getwd(), "/GFF3_additional"))
horseGFF <- file.path(paste0(gffDir, "/ref_EquCab2.0_top_level.gff3"))

# Load previously saved data
load("Globin-Additional_analysis.RData")

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
library(GenomicFeatures)
library(tximport)
library(rjson)
library(reshape2)


# Uncomment functions below to install packages in case you don't have them

# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer", type = "source")
#biocLite("tximport")
#biocLite("GenomicFeatures", type = "source")


# CRAN packages
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("rjson")
#install.packages("reshape2")

###############################################################
# 03 Create NCBI RefSeq TxDb object for transcript annotation #                     #
###############################################################

# Horse NCBI RefSeq .gff3 file used in this analysis
# ftp://ftp.ncbi.nih.gov/genomes/Equus_caballus/GFF/ref_EquCab2.0_top_level.gff3.gz
# ANNOTATION RELEASE NAME: NCBI Equus caballus Annotation Release 102
# ANNOTATION EVIDENCE FREEZE DATE: 19 November 2015
# ANNOTATION RELEASE DATE: 20 November 2015
# ANNOTATION REPORT: http://www.ncbi.nlm.nih.gov/genome/annotation_euk/Equus_caballus/102/
    
# Create DB using .gff3 file
horseDB <- makeTxDbFromGFF(horseGFF, format = "gff3")

# Check databases info
horseDB

################################################################
# 04 Create dataframe of transcripts and gene IDs for tximport #
################################################################

# Check DB key type
keytypes(horseDB)

# Extract keys from DBs
horse_keys <- keys(horseDB, keytype = "GENEID")

# Create transcript-to-gene dataframe for each species
horse_tx2gene <- AnnotationDbi::select(horseDB,
                                       keys    = horse_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

# tximport requires the order to be: transcript name followed by gene ID
horse_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(horse_tx2gene)

##############################################################
# 05 Import salmon TPM estimates and summarise at gene-level #
##############################################################

# Get path to salmon files
horse_files <- list.files(TPMdir, pattern = "quant.sf", full.names = TRUE)
names(horse_files) <- list.files(TPMdir)

# Import transcript-level estimates summarized to the gene-level
# (it's the default when txOut = FALSE or not provided)
# (ignoreTxVersion = TRUE otherwise the horse transcripts names in the UCSC TxDb 
# won't match the ones in salmon's quant.sf)
horse_txi <- tximport(horse_files,
                      type = "salmon",
                      tx2gene = horse_tx2gene)
names(horse_txi)
head(horse_txi$abundance)

############################################
# 06 Remove zero and lowly expressed genes #
############################################

# Convert gene-level TPM abundances into data frames 
horse_TPM <- as.data.frame(horse_txi$abundance) 

# Remove non-expressed genes
horse_nozeros <- horse_TPM[rowSums(horse_TPM) > 0, ]
dim(horse_nozeros)
dim(horse_TPM)

# Remove lowly expressed genes (< 1 TPM in approximately half of total samples)
horse_filt <- horse_nozeros[rowSums(horse_nozeros >= 1) >= 18, ]
dim(horse_filt)
dim(horse_nozeros)

####################################
# 07 Tidy filtered TPM data frames #
####################################

horse_filt %<>% 
    rownames_to_column(var = "Gene_RefSeqID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable) %>% 
    dplyr::mutate(species = "Equine")

# Check reformatted data frames
head(horse_filt)

#######################################
# 08 Convert species column to factor #
#######################################

horse_filt$species %<>%
    factor() 

#####################################################
# 09 Add plotting labels to filtered TPM data frame #
#####################################################

# Create new column
horse_filt$labels <- horse_filt$sample

# Correct plotting labels
horse_filt$labels %<>%
    str_replace("_quant.sf", "") %>% 
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
unique(horse_filt$labels)

# Convert labels into factors and order them
horse_filt$labels %<>%
    factor() %>% 
    forcats::fct_inorder()

# Check factor order
levels(horse_filt$labels)

####################################################
# 10 Add treatment info to filtered TPM data frame #
####################################################

# Create new column
horse_filt$treatment <- horse_filt$labels

# Correct info and convert into factor
horse_filt$treatment %<>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_(D|U)", "Undepleted") %>% 
    factor()

# Check treatment factors
levels(horse_filt$treatment)

#######################
# 11 Save .RData file #
#######################

save.image(file = "Globin-Additional_analysis.RData")

#########################
# 12 Get R session info #
#########################

devtools::session_info()

######################################
# Proceed to Part B of main analysis #
######################################

# File: 07-Globin-Main_analysis-Part_B.R




















