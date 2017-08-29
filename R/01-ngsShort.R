###############################################
# Globin-derived transcripts RNA-seq Analysis #
#       ngsShoRT filtering statistics         #
###############################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: August 25th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Define variables for specific directories
ngsShortDirPE <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/ngsShort/PE"
ngsShortDirSE <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/ngsShort/SE"
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables"

############################################
# 02 Load and/or install required packages #
############################################

# Load packages
library(plyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(devtools)

# Uncomment functions below to install packages in case you don't have them
#install.packages("plyr")
#install.packages("tidyverse")

####################################
# 03 Import ngsShoRT summary files #
####################################

# Create a vector with full paths
summ_filesPE <- list.files(ngsShortDirPE, full.names = TRUE)
summ_filesSE <- list.files(ngsShortDirSE, full.names = TRUE)

# Import paired-end files into a tibble
filtering_stats_PE <- map_df(summ_filesPE, ~ read_table2(.x))
View(filtering_stats_PE)

# Import single-end files into a tibble
filtering_stats_SE <- map_df(summ_filesSE, ~ read_table2(.x))
View(filtering_stats_SE)

########################
# 04 Clean data frames #
########################

# Remove unnecessary columns
filtering_stats_PE %<>%
    dplyr::select(-contains("X"))

filtering_stats_SE %<>%
    dplyr::select(-contains("X"))

# Remove extra symbols and convert data to numeric
filtering_stats_PE$Percent_removed %<>% 
    str_replace("\\(", "") %>% 
    str_replace("\\)", "") %>% 
    str_replace("%", "") %>% 
    as.numeric()

filtering_stats_SE$Percent_removed %<>% 
    str_replace("\\(", "") %>% 
    str_replace("\\)", "") %>% 
    str_replace("%", "") %>% 
    as.numeric()

# Check data frames
View(filtering_stats_PE)
View(filtering_stats_SE)

#############################
# 05 Add labels column (PE) #
#############################

# Add new column for plotting labels
filtering_stats_PE$labels <- filtering_stats_PE$Sample_name

# Correct plotting labels
filtering_stats_PE$labels %<>%
    str_replace("SRR1060762", "Hsa_P10_D") %>% 
    str_replace("SRR1060774", "Hsa_P10_D") %>% 
    str_replace("SRR1060763", "Hsa_P11_D") %>% 
    str_replace("SRR1060775", "Hsa_P11_D") %>%
    str_replace("SRR1060764", "Hsa_P12_D") %>% 
    str_replace("SRR1060776", "Hsa_S01_U") %>% 
    str_replace("SRR1060788", "Hsa_S01_U") %>% 
    str_replace("SRR1060777", "Hsa_S02_U") %>% 
    str_replace("SRR1060789", "Hsa_S02_U") %>% 
    str_replace("SRR1060778", "Hsa_S03_U") %>% 
    str_replace("SRR1060790", "Hsa_S03_U") %>% 
    str_replace("SRR1060779", "Hsa_S04_U") %>% 
    str_replace("SRR1060791", "Hsa_S04_U") %>% 
    str_replace("SRR1060780", "Hsa_S05_U") %>% 
    str_replace("SRR1060792", "Hsa_S05_U") %>% 
    str_replace("SRR1060781", "Hsa_S06_U") %>% 
    str_replace("SRR1060793", "Hsa_S06_U") %>% 
    str_replace("SRR1060782", "Hsa_P19_U") %>% 
    str_replace("SRR1060794", "Hsa_P19_U") %>% 
    str_replace("SRR1060783", "Hsa_P20_U") %>% 
    str_replace("SRR1060795", "Hsa_P20_U") %>% 
    str_replace("SRR1060784", "Hsa_P21_U") %>% 
    str_replace("SRR1060796", "Hsa_P21_U") %>% 
    str_replace("SRR1060785", "Hsa_P22_U") %>% 
    str_replace("SRR1060797", "Hsa_P22_U") %>%
    str_replace("SRR1060786", "Hsa_P23_U") %>% 
    str_replace("SRR1060798", "Hsa_P23_U") %>% 
    str_replace("SRR1060787", "Hsa_P24_U") %>% 
    str_replace("SRR1060799", "Hsa_P24_U") %>%
    str_replace("SRR1060753", "Hsa_S01_D") %>% 
    str_replace("SRR1060765", "Hsa_S01_D") %>% 
    str_replace("SRR1060754", "Hsa_S02_D") %>% 
    str_replace("SRR1060766", "Hsa_S02_D") %>% 
    str_replace("SRR1060755", "Hsa_S03_D") %>% 
    str_replace("SRR1060767", "Hsa_S03_D") %>%
    str_replace("SRR1060756", "Hsa_S04_D") %>% 
    str_replace("SRR1060768", "Hsa_S04_D") %>% 
    str_replace("SRR1060757", "Hsa_S05_D") %>% 
    str_replace("SRR1060769", "Hsa_S05_D") %>% 
    str_replace("SRR1060758", "Hsa_S06_D") %>% 
    str_replace("SRR1060770", "Hsa_S06_D") %>% 
    str_replace("SRR1060759", "Hsa_P07_D") %>% 
    str_replace("SRR1060771", "Hsa_P07_D") %>%
    str_replace("SRR1060760", "Hsa_P08_D") %>% 
    str_replace("SRR1060772", "Hsa_P08_D") %>%
    str_replace("SRR1060761", "Hsa_P09_D") %>% 
    str_replace("SRR1060773", "Hsa_P09_D") %>%
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
    str_replace("_W-1_F_00\\d", "") %>%
    str_replace("6511", "Bta_01_U") %>%
    str_replace("6514", "Bta_02_U") %>%
    str_replace("6520", "Bta_03_U") %>%
    str_replace("6522", "Bta_04_U") %>%
    str_replace("6526", "Bta_05_U") %>%
    str_replace("6635", "Bta_06_U") %>%
    str_replace("6636", "Bta_07_U") %>%
    str_replace("6637", "Bta_08_U") %>%
    str_replace("6644", "Bta_09_U") %>%
    str_replace("6698", "Bta_10_U")

# Convert labels into factors and order them
filtering_stats_PE$labels %<>%
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
                      "Bta_01_U", "Bta_02_U", "Bta_03_U", "Bta_04_U",
                      "Bta_05_U", "Bta_06_U", "Bta_07_U", "Bta_08_U",
                      "Bta_09_U", "Bta_10_U"))

# Check factor order
levels(filtering_stats_PE$labels)

# Check data frame
View(filtering_stats_PE)

#############################
# 06 Add labels column (SE) #
#############################

# Add new column for plotting labels
filtering_stats_SE$labels <- filtering_stats_SE$Sample_name

# Correct plotting labels
filtering_stats_SE$labels %<>%
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

# Convert labels into factors and order them
filtering_stats_SE$labels %<>%
    factor(levels = c("Eca_T01_U", "Eca_T02_U", "Eca_T03_U", "Eca_T04_U",
                      "Eca_T05_U", "Eca_T06_U", "Eca_T07_U", "Eca_T08_U",
                      "Eca_T09_U", "Eca_T10_U", "Eca_T11_U", "Eca_T12_U",
                      "Eca_S13_U", "Eca_T14_U", "Eca_T15_U", "Eca_T16_U",
                      "Eca_T17_U", "Eca_T18_U", "Eca_T19_U", "Eca_T20_U",
                      "Eca_T21_U", "Eca_T22_U", "Eca_T23_U", "Eca_T24_U",
                      "Eca_T25_U", "Eca_T26_U", "Eca_T27_U", "Eca_T28_U",
                      "Eca_T29_U", "Eca_T30_U", "Eca_T31_U", "Eca_T32_U",
                      "Eca_S33_U", "Eca_S34_U", "Eca_S35_U", "Eca_S36_U",
                      "Eca_S37_U"))

# Check factor order
levels(filtering_stats_SE$labels)

# Check data frame
View(filtering_stats_SE)

###########################
# 07 Add treatment column #
###########################

# Add new column for traetment
filtering_stats_PE$treatment <- filtering_stats_PE$labels
filtering_stats_SE$treatment <- c("Undepleted")

# Correct values
filtering_stats_PE$treatment %<>% 
    stringr::str_replace("Bta_\\d\\d_", "") %>% 
    stringr::str_replace("Hsa_(S|P)\\d\\d_", "") %>% 
    stringr::str_replace("Ssc_\\d\\d_", "") %>%
    stringr::str_replace("U", "Undepleted") %>% 
    stringr::str_replace("D", "Globin depleted")

# Check data frames
View(filtering_stats_PE)
View(filtering_stats_SE)

#########################
# 08 Add species column #
#########################

# Add new column for species
filtering_stats_PE$species <- filtering_stats_PE$labels
filtering_stats_SE$species <- c("Equine")

# Correct values
filtering_stats_PE$species %<>% 
    stringr::str_replace("Bta_\\d\\d_(U|D)", "Bovine") %>% 
    stringr::str_replace("Hsa_(S|P)\\d\\d_(U|D)", "Human") %>% 
    stringr::str_replace("Ssc_\\d\\d_(U|D)", "Porcine")

# Check data frames
View(filtering_stats_PE)
View(filtering_stats_SE)

########################################
# 09 Order data frames and export data #
########################################

# Order columns
filtering_stats_PE %<>% 
    dplyr::select(species, labels, treatment, everything())

filtering_stats_SE %<>% 
    dplyr::select(species, labels, treatment, everything())

# Order rows by labels
filtering_stats_PE %<>% 
    dplyr::arrange(labels)

filtering_stats_SE %<>% 
    dplyr::arrange(labels)

# Export tidy summary stats
write_csv(filtering_stats_PE,
          path = file.path(paste0(tablesDir, "/ngsShort_stats_PE.csv")),
          col_names = TRUE)

write_csv(filtering_stats_SE,
          path = file.path(paste0(tablesDir, "/ngsShort_stats_SE.csv")),
          col_names = TRUE)

############################
# 08 R session information #
############################

devtools::session_info()


