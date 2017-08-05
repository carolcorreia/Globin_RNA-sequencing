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
cattle_nozeros <- cattle_TPM[rowSums(cattle_TPM) > 0, ]
dim(cattle_nozeros)
dim(cattle_TPM)

horse_nozeros <- horse_TPM[rowSums(horse_TPM) > 0, ]
dim(horse_nozeros)
dim(horse_TPM)

human_nozeros <- human_TPM[rowSums(human_TPM) > 0, ]
dim(human_nozeros)
dim(human_TPM)

pig_nozeros <- pig_TPM[rowSums(pig_TPM) > 0, ]
dim(pig_nozeros)
dim(pig_TPM)

# Remove lowly expressed genes (< 1 TPM in one treatment group)
cattle_filt <- cattle_nozeros[rowSums(cattle_nozeros >= 1) >= 10, ]
dim(cattle_filt)
dim(cattle_nozeros)

horse_filt <- horse_nozeros[rowSums(horse_nozeros >= 1) >= 37, ]
dim(horse_filt)
dim(horse_nozeros)

human_filt <- human_nozeros[rowSums(human_nozeros >= 1) >= 12, ]
dim(human_filt)
dim(human_nozeros)

pig_filt <- pig_nozeros[rowSums(pig_nozeros >= 1) >= 12, ]
dim(pig_filt)
dim(pig_nozeros)

####################################
# 08 Tidy filtered TPM data frames #
####################################

cattle_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

horse_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

human_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

pig_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

# Check reformatted data frames
head(cattle_filt)
head(horse_filt)
head(human_filt)
head(pig_filt)

#########################################################
# 09 Row-bind filtered TPM data frames from all species #
#########################################################

TPM_filt_all <- dplyr::bind_rows(cattle_filt, horse_filt)
TPM_filt_all <- dplyr::bind_rows(TPM_filt_all, human_TPM_filt)

cattle_filt %>% 
    dplyr::bind_rows(horse_filt) %>% 
    dplyr::bind_rows(human_filt) %>% 
    dplyr::bind_rows(pig_filt) %>% 
    as.tibble() -> TPM_filt_all

# Check total number of rows
dim(TPM_filt_all)
dim(cattle_filt) + dim(horse_filt) + dim(human_filt) + dim(pig_filt)

#####################################################
# 10 Add plotting labels to filtered TPM data frame #
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
    str_replace("SRR3671009", "Eca_T01") %>%
    str_replace("SRR3671010", "Eca_T02") %>%
    str_replace("SRR3671011", "Eca_T03") %>%
    str_replace("SRR3671012", "Eca_T04") %>%
    str_replace("SRR3671013", "Eca_T05") %>%
    str_replace("SRR3671014", "Eca_T06") %>%
    str_replace("SRR3671015", "Eca_T07") %>%
    str_replace("SRR3671016", "Eca_T08") %>%
    str_replace("SRR3671017", "Eca_T09") %>%
    str_replace("SRR3671018", "Eca_T10") %>%
    str_replace("SRR3671019", "Eca_T11") %>%
    str_replace("SRR3671020", "Eca_T12") %>%
    str_replace("SRR3671021", "Eca_S13") %>%
    str_replace("SRR3671022", "Eca_T14") %>%
    str_replace("SRR3671023", "Eca_T15") %>%
    str_replace("SRR3671024", "Eca_T16") %>%
    str_replace("SRR3671025", "Eca_T17") %>%
    str_replace("SRR3671026", "Eca_T18") %>%
    str_replace("SRR3671027", "Eca_T19") %>%
    str_replace("SRR3671028", "Eca_T20") %>%
    str_replace("SRR3671029", "Eca_T21") %>%
    str_replace("SRR3671030", "Eca_T22") %>%
    str_replace("SRR3671031", "Eca_T23") %>%
    str_replace("SRR3671032", "Eca_T24") %>%
    str_replace("SRR3671033", "Eca_T25") %>%
    str_replace("SRR3671034", "Eca_T26") %>%
    str_replace("SRR3671035", "Eca_T27") %>%
    str_replace("SRR3671036", "Eca_T28") %>%
    str_replace("SRR3671037", "Eca_T29") %>%
    str_replace("SRR3671038", "Eca_T30") %>%
    str_replace("SRR3671039", "Eca_T31") %>%
    str_replace("SRR3671040", "Eca_T32") %>%
    str_replace("SRR3671041", "Eca_S33") %>%
    str_replace("SRR3671042", "Eca_S34") %>%
    str_replace("SRR3671043", "Eca_S35") %>%
    str_replace("SRR3671044", "Eca_S36") %>%
    str_replace("SRR3671045", "Eca_S37")

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
# 11 Add treatment info to filtered TPM data frame #
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

##################################################
# 12 Add species info to filtered TPM data frame #
##################################################

# Add species column A_TPM_filt
TPM_filt_all$species <- TPM_filt_all$Ensembl_gene_ID

# Correct info and convert into factor
TPM_filt_all$species %<>%
    stringr::str_replace("ENSG0.+", "Human") %>% 
    stringr::str_replace("ENSSSCG0.+", "Porcine") %>% 
    stringr::str_replace("ENSECAG0.+", "Equine") %>%
    stringr::str_replace("ENSBTAG0.+", "Bovine") %>% 
    factor(levels = c("Human", "Porcine", "Equine", "Bovine"))

# Check species factors
levels(TPM_filt_all$species)

#####################
#  Save .RData file #
#####################

save.image(file = "Globin-RNA-seqAnalysis.RData")

#######################
#  Get R session info #
#######################

devtools::session_info()
























