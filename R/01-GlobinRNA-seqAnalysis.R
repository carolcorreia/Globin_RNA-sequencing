# HBA and HBB transcripts RNA-seq Analysis
# Author: Carolina C. Correia
# Date: August 4th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")
getwd()

# Define variables for specific directories
workDir <- getwd()
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"
gtfDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R/GTF_ensembl"
cattleDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/TPM/cattle_TPM"
horseDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/TPM/horse_TPM"
humanDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/TPM/human_TPM"
pigDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/TPM/pig_TPM"

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

# Check DB key type
keytypes(cattleDB)
keytypes(horseDB)
keytypes(humanDB)
keytypes(pigDB)

# Extract keys from DBs
cattle_keys <- biomaRt::keys(cattleDB, keytype = "GENEID")
horse_keys <- biomaRt::keys(horseDB, keytype = "GENEID")
human_keys <- biomaRt::keys(humanDB, keytype = "GENEID")
pig_keys <- biomaRt::keys(pigDB, keytype = "GENEID")

# Create transcript-to-gene df for each species
cattle_tx2gene <- AnnotationDbi::select(cattleDB,
                                        keys = cattle_keys,
                                        keytype = "GENEID",
                                        columns = "TXNAME")

horse_tx2gene <- AnnotationDbi::select(horseDB,
                                       keys = horse_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

human_tx2gene <- AnnotationDbi::select(humanDB,
                                       keys = human_keys,
                                       keytype = "GENEID",
                                       columns = "TXNAME")

pig_tx2gene <- AnnotationDbi::select(pigDB,
                                     keys = pig_keys,
                                     keytype = "GENEID",
                                     columns = "TXNAME")

# tximport requires the order to be: transcript name followed by gene ID
cattle_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(cattle_tx2gene)


horse_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(horse_tx2gene)


human_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(human_tx2gene)


pig_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)
head(pig_tx2gene)

##############################################################
# 05 Import salmon TPM estimates and summarise at gene-level #
##############################################################

# Get paths to salmon files
cattle_files <- list.files(cattleDir, pattern = "quant.sf", full.names = TRUE)
names(cattle_files) <- list.files(cattleDir)

horse_files <- list.files(horseDir, pattern = "quant.sf", full.names = TRUE)
names(horse_files) <- list.files(horseDir)

human_files <- list.files(humanDir, pattern = "quant.sf", full.names = TRUE)
names(human_files) <- list.files(humanDir)

pig_files <- list.files(pigDir, pattern = "quant.sf", full.names = TRUE)
names(pig_files) <- list.files(pigDir)

# Import transcript-level estimates summarized to the gene-level
# (it's the default when txOut = FALSE or not provided)
# (ignoreTxVersion = TRUE otherwise the transcripts names in the GTF file 
# won't match the ones in salmon's quant.sf)
cattle_txi <- tximport(cattle_files,
                       type = "salmon",
                       tx2gene = cattle_tx2gene,
                       ignoreTxVersion = TRUE)
names(cattle_txi)
head(cattle_txi$abundance)


horse_txi <- tximport(horse_files,
                      type = "salmon",
                      tx2gene = horse_tx2gene,
                      ignoreTxVersion = TRUE)
names(horse_txi)
head(horse_txi$abundance)


human_txi <- tximport(human_files,
                      type = "salmon",
                      tx2gene = human_tx2gene,
                      ignoreTxVersion = TRUE)
names(human_txi)
head(human_txi$abundance)


pig_txi <- tximport(pig_files,
                    type = "salmon",
                    tx2gene = pig_tx2gene,
                    ignoreTxVersion = TRUE)
names(pig_txi)
head(pig_txi$abundance)

#####################
# 06 Get gene names #
#####################

# Get datasets names from Ensembl biomaRt
listMarts(mart = NULL, host="www.ensembl.org", path="/biomart/martservice",
          port = 80)
mart_datasets <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart_datasets)

# Convert gene-level TPM abundances into data frames
cattle_TPM <- as.data.frame(cattle_txi$abundance)
horse_TPM <- as.data.frame(horse_txi$abundance)
human_TPM <- as.data.frame(human_txi$abundance)
pig_TPM <- as.data.frame(pig_txi$abundance)

# Get gene symbols for cattle
mart_cattle = useMart("ENSEMBL_MART_ENSEMBL",
                      dataset = "btaurus_gene_ensembl")
listFilters(mart_cattle)

cattle_genes <- getBM(attributes = c("ensembl_gene_id",
                                     "external_gene_name"),
                     values = rownames(cattle_TPM),
                     mart = mart_cattle)
head(cattle_genes)

# Get gene symbols for horse
mart_horse = useMart("ENSEMBL_MART_ENSEMBL",
                     dataset = "ecaballus_gene_ensembl")
listFilters(mart_horse)

horse_genes <- getBM(attributes = c("ensembl_gene_id",
                                    "external_gene_name"),
                     values = rownames(horse_TPM),
                     mart = mart_horse)
head(horse_genes)

# Get gene symbols for human
mart_human = useMart("ENSEMBL_MART_ENSEMBL",
                     dataset = "hsapiens_gene_ensembl")
listFilters(mart_human)

human_genes <- getBM(attributes = c("ensembl_gene_id",
                                    "external_gene_name"),
                     values = rownames(human_TPM),
                     mart = mart_human)
head(human_genes)

# Get gene symbols for pig
mart_pig = useMart("ENSEMBL_MART_ENSEMBL",
                   dataset = "sscrofa_gene_ensembl")
listFilters(mart_pig)

pig_genes <- getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   values = rownames(pig_TPM),
                   mart = mart_pig)
head(pig_genes)

############################################
# 07 Remove zero and lowly expressed genes #
############################################

# Remove non-expressed genes
h_TPM_nozeros <- human_TPM[rowSums(human_TPM) > 0, ]
dim(h_TPM_nozeros)
dim(human_TPM)

p_TPM_nozeros <- pig_TPM[rowSums(pig_TPM) > 0, ]
dim(p_TPM_nozeros)
dim(pig_TPM)

c_TPM_nozeros <- cattle_TPM[rowSums(cattle_TPM) > 0, ]
dim(c_TPM_nozeros)
dim(cattle_TPM)

# Remove lowly expressed genes (< 1 TPM)
h_TPM_filt <- h_TPM_nozeros[rowSums(h_TPM_nozeros >= 1) >= 12, ]
dim(h_TPM_filt)
dim(h_TPM_nozeros)

p_TPM_filt <- p_TPM_nozeros[rowSums(p_TPM_nozeros >= 1) >= 12, ]
dim(p_TPM_filt)
dim(p_TPM_nozeros)

c_TPM_filt <- c_TPM_nozeros[rowSums(c_TPM_nozeros >= 1) >= 10, ]
dim(c_TPM_filt)
dim(c_TPM_nozeros)


# Tidy filt TPM dfs
h_TPM_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(h_TPM_filt)

p_TPM_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(p_TPM_filt)

c_TPM_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(c_TPM_filt)


#####################
#  Save .RData file #
#####################

save.image(file = "Globin-RNA-seqAnalysis.RData")

#######################
#  Get R session info #
#######################

devtools::session_info()
























