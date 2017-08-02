# HBA and HBB transcripts RNA-seq Analysis

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")
getwd()

# Define variables for working and image directories
workDir <- getwd()
img_dir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"

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

############################################
# 03  #
############################################





























