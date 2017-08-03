# HBA and HBB transcripts RNA-seq Analysis

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")
getwd()

# Define variables for specific directories
workDir <- getwd()
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"
tpmDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/TPM"
gtfDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R/GTF_ensembl"

# Load previously saved data
load("Globin-RNA-seqAnalysis.RData")

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(plyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(repurrrsive)
library(ggrepel)
library(ggpubr)
library(ggjoy)
library(waffle)
library(rtracklayer)
library(GenomicFeatures)
library(tximport)
library(rjson)
library(reshape2)
library(biomaRt)
library(viridis)

# Uncomment functions below to install packages in case you don't have them

# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer", type = "source")
#biocLite("tximport")
#biocLite("GenomicFeatures", type = "source")
#biocLite("biomaRt")


# CRAN packages
#install.packages("devtools")
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("rjson")
#install.packages("ggrepel")
#install.packages("ggpubr")
#install.packages("ggjoy")
#install.packages("reshape2")
#install.packages("viridis")
#devtools::install_github("hrbrmstr/waffle")
#devtools::install_github("jennybc/repurrrsive")

##################################################
# 03 Create annotation DB from Ensembl GTF files #
##################################################

# Ensembl Release 88 GTF files were downloaded from the URLs below:
# ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz
# ftp://ftp.ensembl.org/pub/release-88/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.88.gtf.gz
# ftp://ftp.ensembl.org/pub/release-88/gtf/bos_taurus/Bos_taurus.UMD3.1.88.gtf.gz
# ftp://ftp.ensembl.org/pub/release-88/gtf/equus_caballus/Equus_caballus.EquCab2.88.gtf.gz


# Get full paths for GTF files
gtf_files <- list.files(gtfDir, pattern = ".gtf", full.names = TRUE)
gtf_files

# Create annotation DBs for each species
cattleDB <- makeTxDbFromGFF(file = gtf_files[1],
                            format = "gtf",
                            dataSource = "Ensembl GTF Bos_taurus.UMD3.1.88",
                            organism = "Bos taurus")

horseDB <- makeTxDbFromGFF(file = gtf_files[2],
                           format = "gtf",
                           dataSource = "Ensembl GTF Equus_caballus.EquCab2.88",
                           organism = "Equus caballus")

humanDB <- makeTxDbFromGFF(file = gtf_files[3],
                           format = "gtf",
                           dataSource = "Ensembl GTF Homo_sapiens.GRCh38.88",
                           organism = "Homo sapiens")

pigDB <- makeTxDbFromGFF(file = gtf_files[4],
                         format = "gtf",
                         dataSource = "Ensembl GTF Sus_scrofa.Sscrofa10.2.88",
                         organism = "Sus scrofa")


# Check DB info
cattleDB
horseDB
humanDB
pigDB


################################################################
# 04 Create dataframe of transcripts and gene IDs for tximport #
################################################################


keytypes(humanDB)
human_keys <- keys(humanDB, keytype = "GENEID")
human_tx2gene <- AnnotationDbi::select(humanDB,
                                       keys = human_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

list_DBs <- list(cattleDB, horseDB, humanDB, pigDB)
keys_DBs <- list(cattle_keys, horse_keys, human_keys, pig_keys)

map_df(list_DBs, AnnotationDbi::select(keys = keys_DBs,
                                        keytype = "GENEID",
                                        columns = "TXNAME"))

# tximport requires the order to be: transcript name followed by gene ID
human_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)

head(human_tx2gene)

horse_keys <- keys(horseDB, keytype = "GENEID")



keytypes(pigDB)
pig_keys <- keys(pigDB, keytype = "GENEID")
P_tx2gene <- AnnotationDbi::select(pigDB,
                                   keys = pig_keys,
                                   keytype = "GENEID",
                                   columns = "TXNAME")

P_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)

head(P_tx2gene)

keytypes(cattleDB)
cattle_keys <- keys(cattleDB, keytype = "GENEID")
C_tx2gene <- AnnotationDbi::select(cattleDB,
                                   keys = cattle_keys,
                                   keytype = "GENEID",
                                   columns = "TXNAME")

C_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)

head(C_tx2gene)

########################
# 05 Summarise gene-level TPM #
########################

# Get directories' path
path_dirs <- list.dirs(tpmDir)

# Get file names
salmon_files <- list.files(path_dirs, pattern = "quant.sf", full.names = TRUE)








save.image(file = "Globin-RNA-seqAnalysis.RData")

devtools::session_info()
























