# HBA1, HBA2, and HBB transcripts RNA-seq Analysis - Part 1
# Author: Carolina C. Correia
# Date: August 16th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Define variables for specific directories
salmonDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon"
humanDir <- file.path(paste0(salmonDir, "/TPM/human_TPM"))
pigDir <- file.path(paste0(salmonDir, "/TPM/pig_TPM"))
horseDir <- file.path(paste0(salmonDir, "/UCSC_TPM/horse_TPM"))
cattleDir <- file.path(paste0(salmonDir, "/TPM/cattle_TPM"))

gffDir <- file.path(paste0(getwd(), "/GFF3_RefSeq"))
humanGFF <- file.path(paste0(gffDir, "/ref_GRCh38.p7_top_level.gff3"))
pigGFF <- file.path(paste0(gffDir, "/ref_Sscrofa11.1_top_level.gff3"))
cattleGFF <- file.path(paste0(gffDir, "/ref_Bos_taurus_UMD_3.1.1_top_level.gff3")) 
    
# Load previously saved data
load("Globin-RNA-seqAnalysis.RData")

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
#biocLite("TxDb.Btaurus.UCSC.bosTau8.refGene")
#biocLite("TxDb.Sscrofa.UCSC.susScr3.refGene")


# CRAN packages
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("rjson")
#install.packages("reshape2")

################################################################
# 03 Create NCBI RefSeq TxDb objects for transcript annotation #
#                      (human, pig, cow)                       #
################################################################

# Human NCBI RefSeq .gff3 file used in this analysis
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz
# ANNOTATION RELEASE NAME: NCBI Homo sapiens Annotation Release 108
# ANNOTATION EVIDENCE FREEZE DATE: 5 May 2016
# ANNOTATION RELEASE DATE: 6 June 2016
# ANNOTATION REPORT: http://www.ncbi.nlm.nih.gov/genome/annotation_euk/Homo_sapiens/108/
    
# Pig NCBI RefSeq .gff3 file used in this analysis
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Sus_scrofa/GFF/ref_Sscrofa11.1_top_level.gff3.gz
# ANNOTATION RELEASE NAME: NCBI Sus scrofa Annotation Release 106
# ANNOTATION EVIDENCE FREEZE DATE: 3 May 2017
# ANNOTATION RELEASE DATE: 13 May 2017
# ANNOTATION REPORT: https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Sus_scrofa/106/

# Cow NCBI RefSeq .gff3 file used in this analysis
# ftp://ftp.ncbi.nih.gov/genomes/Bos_taurus/GFF/ref_Bos_taurus_UMD_3.1.1_top_level.gff3.gz
# ANNOTATION RELEASE NAME: NCBI Bos taurus Annotation Release 105
# ANNOTATION EVIDENCE FREEZE DATE: 20 January 2016
# ANNOTATION RELEASE DATE: 26 January 2016
# ANNOTATION REPORT: http://www.ncbi.nlm.nih.gov/genome/annotation_euk/Bos_taurus/105/
    
# Create DBs using .gff3 files
humanDB <- makeTxDbFromGFF(humanGFF, format = "gff3")
pigDB <- makeTxDbFromGFF(pigGFF, format = "gff3")
cattleDB <- makeTxDbFromGFF(cattleGFF, format = "gff3")

# Check databases info
humanDB
pigDB
cattleDB

################################################################
# 04 Create UCSC TxDb object for transcript annotation (horse) #
################################################################

# UCSC genome refGene (NCBI RefSeq genes track) table used: 
# Sep. 2007 (equCab2, Broad Institute EquCab2) assembly of the horse genome

# Display the list of tables known to work with makeTxDbFromUCSC()
supportedUCSCtables()

# Retrieve a full transcript dataset from UCSC
horseDB <- makeTxDbFromUCSC(genome = "equCab2", tablename = "refGene")

# Check DB info
horseDB

################################################################
# 05 Create dataframe of transcripts and gene IDs for tximport #
################################################################

# Check DB key type
keytypes(humanDB)
keytypes(pigDB)
keytypes(horseDB)
keytypes(cattleDB)

# Extract keys from DBs
human_keys <- keys(humanDB, keytype = "GENEID")
pig_keys <- keys(pigDB, keytype = "GENEID")
horse_keys <- keys(horseDB, keytype = "GENEID")
cattle_keys <- keys(cattleDB, keytype = "GENEID")

# Create transcript-to-gene dataframe for each species
human_tx2gene <- AnnotationDbi::select(humanDB,
                                       keys = human_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

pig_tx2gene <- AnnotationDbi::select(pigDB,
                                     keys = pig_keys,
                                     keytype = "GENEID",
                                     columns = "TXNAME")

horse_tx2gene <- AnnotationDbi::select(horseDB,
                                       keys = horse_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

cattle_tx2gene <- AnnotationDbi::select(cattleDB,
                                        keys = cattle_keys,
                                        keytype = "GENEID",
                                        columns = "TXNAME")

# tximport requires the order to be: transcript name followed by gene ID
human_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(human_tx2gene)


pig_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(pig_tx2gene)


horse_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(horse_tx2gene)


cattle_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(cattle_tx2gene)

##############################################################
# 06 Import salmon TPM estimates and summarise at gene-level #
##############################################################

# Get paths to salmon files
human_files <- list.files(humanDir, pattern = "quant.sf", full.names = TRUE)
names(human_files) <- list.files(humanDir)

pig_files <- list.files(pigDir, pattern = "quant.sf", full.names = TRUE)
names(pig_files) <- list.files(pigDir)

horse_files <- list.files(horseDir, pattern = "quant.sf", full.names = TRUE)
names(horse_files) <- list.files(horseDir)

cattle_files <- list.files(cattleDir, pattern = "quant.sf", full.names = TRUE)
names(cattle_files) <- list.files(cattleDir)

# Import transcript-level estimates summarized to the gene-level
# (it's the default when txOut = FALSE or not provided)
# (ignoreTxVersion = TRUE otherwise the horse transcripts names in the UCSC TxDb 
# won't match the ones in salmon's quant.sf)
human_txi <- tximport(human_files,
                      type = "salmon",
                      tx2gene = human_tx2gene)
names(human_txi)
head(human_txi$abundance)


pig_txi <- tximport(pig_files,
                    type = "salmon",
                    tx2gene = pig_tx2gene)
names(pig_txi)
head(pig_txi$abundance)


horse_txi <- tximport(horse_files,
                      type = "salmon",
                      tx2gene = horse_tx2gene,
                      ignoreTxVersion = TRUE)
names(horse_txi)
head(horse_txi$abundance)


cattle_txi <- tximport(cattle_files,
                       type = "salmon",
                       tx2gene = cattle_tx2gene)
names(cattle_txi)
head(cattle_txi$abundance)

############################################
# 07 Remove zero and lowly expressed genes #
############################################

# Convert gene-level TPM abundances into data frames 
human_TPM <- as.data.frame(human_txi$abundance) 
pig_TPM <- as.data.frame(pig_txi$abundance)
horse_TPM <- as.data.frame(horse_txi$abundance)
cattle_TPM <- as.data.frame(cattle_txi$abundance)

# Remove non-expressed genes
human_nozeros <- human_TPM[rowSums(human_TPM) > 0, ]
dim(human_nozeros)
dim(human_TPM)

pig_nozeros <- pig_TPM[rowSums(pig_TPM) > 0, ]
dim(pig_nozeros)
dim(pig_TPM)

horse_nozeros <- horse_TPM[rowSums(horse_TPM) > 0, ]
dim(horse_nozeros)
dim(horse_TPM)

cattle_nozeros <- cattle_TPM[rowSums(cattle_TPM) > 0, ]
dim(cattle_nozeros)
dim(cattle_TPM)

# Remove lowly expressed genes (< 1 TPM in one treatment group)
human_filt <- human_nozeros[rowSums(human_nozeros >= 1) >= 12, ]
dim(human_filt)
dim(human_nozeros)

pig_filt <- pig_nozeros[rowSums(pig_nozeros >= 1) >= 12, ]
dim(pig_filt)
dim(pig_nozeros)

horse_filt <- horse_nozeros[rowSums(horse_nozeros >= 1) >= 37, ]
dim(horse_filt)
dim(horse_nozeros)

cattle_filt <- cattle_nozeros[rowSums(cattle_nozeros >= 1) >= 10, ]
dim(cattle_filt)
dim(cattle_nozeros)

####################################
# 08 Tidy filtered TPM data frames #
####################################
human_filt %<>% 
    rownames_to_column(var = "Gene_RefSeqID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable) %>% 
    dplyr::mutate(species = "Human")

pig_filt %<>% 
    rownames_to_column(var = "Gene_RefSeqID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable) %>% 
    dplyr::mutate(species = "Porcine")

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
head(human_filt)
head(pig_filt)
head(horse_filt)
head(cattle_filt)

#########################################################
# 09 Row-bind filtered TPM data frames from all species #
#########################################################

human_filt %>% 
dplyr::bind_rows(pig_filt) %>% 
dplyr::bind_rows(horse_filt) %>% 
dplyr::bind_rows(cattle_filt) %>% 
    as.tibble() -> TPM_filt_all

# Check total number of rows
length(TPM_filt_all$Gene_RefSeqID) ==
    length(human_filt$Gene_RefSeqID) +
    length(pig_filt$Gene_RefSeqID) +
    length(horse_filt$Gene_RefSeqID) +
    length(cattle_filt$Gene_RefSeqID)

# Visualise TPM df
TPM_filt_all

#######################################
# 10 Convert species column to factor #
#######################################

# Order of the levels match the order of the first appearance in the data
TPM_filt_all$species %<>%
    factor() %>%
    forcats::fct_inorder()

# Check species factors
levels(TPM_filt_all$species)

#####################################################
# 11 Add plotting labels to filtered TPM data frame #
#####################################################

# Create new column
TPM_filt_all$labels <- TPM_filt_all$sample

# Correct plotting labels
TPM_filt_all$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("NGD_(FALSE|TRUE)_", "") %>% 
    str_replace("GD_(FALSE|TRUE)_", "") %>% 
    str_replace("Subj10", "Hsa_P10_D") %>% 
    str_replace("Subj11", "Hsa_P11_D") %>% 
    str_replace("Subj12", "Hsa_P12_D") %>% 
    str_replace("Subj13", "Hsa_S01_U") %>% 
    str_replace("Subj14", "Hsa_S02_U") %>% 
    str_replace("Subj15", "Hsa_S03_U") %>% 
    str_replace("Subj16", "Hsa_S04_U") %>% 
    str_replace("Subj17", "Hsa_S05_U") %>% 
    str_replace("Subj18", "Hsa_S06_U") %>% 
    str_replace("Subj19", "Hsa_P19_U") %>% 
    str_replace("Subj20", "Hsa_P20_U") %>% 
    str_replace("Subj21", "Hsa_P21_U") %>% 
    str_replace("Subj22", "Hsa_P22_U") %>% 
    str_replace("Subj23", "Hsa_P23_U") %>% 
    str_replace("Subj24", "Hsa_P24_U") %>% 
    str_replace("Subj1", "Hsa_S01_D") %>% 
    str_replace("Subj2", "Hsa_S02_D") %>% 
    str_replace("Subj3", "Hsa_S03_D") %>% 
    str_replace("Subj4", "Hsa_S04_D") %>% 
    str_replace("Subj5", "Hsa_S05_D") %>% 
    str_replace("Subj6", "Hsa_S06_D") %>% 
    str_replace("Subj7", "Hsa_P07_D") %>% 
    str_replace("Subj8", "Hsa_P08_D") %>% 
    str_replace("Subj9", "Hsa_P09_D") %>% 
    str_replace("7197", "Ssc_01") %>% 
    str_replace("7199", "Ssc_02") %>% 
    str_replace("7210", "Ssc_03") %>% 
    str_replace("7312", "Ssc_04") %>% 
    str_replace("7349", "Ssc_05") %>%
    str_replace("7413", "Ssc_06") %>%
    str_replace("7437", "Ssc_07") %>%
    str_replace("7439", "Ssc_08") %>%
    str_replace("7467", "Ssc_09") %>%
    str_replace("7468", "Ssc_10") %>%
    str_replace("7472", "Ssc_11") %>%
    str_replace("7474", "Ssc_12") %>%
    str_replace("(CT|C)", "_U") %>% 
    str_replace("(GD|T)", "_D") %>% 
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
    factor(levels = c("Hsa_S01_U", "Hsa_S02_U", "Hsa_S03_U", "Hsa_S04_U",
                      "Hsa_S05_U", "Hsa_S06_U", "Hsa_P19_U", "Hsa_P20_U",
                      "Hsa_P21_U", "Hsa_P22_U", "Hsa_P23_U", "Hsa_P24_U",
                      "Hsa_S01_D", "Hsa_S02_D", "Hsa_S03_D", "Hsa_S04_D",
                      "Hsa_S05_D", "Hsa_S06_D", "Hsa_P07_D", "Hsa_P08_D",
                      "Hsa_P09_D", "Hsa_P10_D", "Hsa_P11_D", "Hsa_P12_D",
                      "Ssc_01_U", "Ssc_02_U", "Ssc_03_U", "Ssc_04_U",
                      "Ssc_05_U", "Ssc_06_U", "Ssc_07_U", "Ssc_08_U",
                      "Ssc_09_U", "Ssc_10_U", "Ssc_11_U", "Ssc_12_U",
                      "Ssc_01_D", "Ssc_02_D", "Ssc_03_D", "Ssc_04_D",
                      "Ssc_05_D", "Ssc_06_D", "Ssc_07_D", "Ssc_08_D",
                      "Ssc_09_D", "Ssc_10_D", "Ssc_11_D", "Ssc_12_D",
                      "Eca_T01_U", "Eca_T02_U", "Eca_T03_U", "Eca_T04_U",
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
# 12 Add treatment info to filtered TPM data frame #
####################################################

# Create new column
TPM_filt_all$treatment <- TPM_filt_all$labels

# Correct info and convert into factor
TPM_filt_all$treatment %<>% 
    stringr::str_replace("Bta_\\d\\d_", "") %>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_", "") %>% 
    stringr::str_replace("Hsa_(S|P)\\d\\d_", "") %>% 
    stringr::str_replace("Ssc_\\d\\d_", "") %>%
    stringr::str_replace("U", "Undepleted") %>% 
    stringr::str_replace("D", "Globin depleted") %>%
    factor(levels = c("Undepleted", "Globin depleted"))

# Check treatment factors
levels(TPM_filt_all$treatment)

#######################
# 13 Save .RData file #
#######################

save.image(file = "Globin-RNA-seqAnalysis.RData")

#########################
# 14 Get R session info #
#########################

devtools::session_info()

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02-GlobinRNA-seqAnalysis.R




















