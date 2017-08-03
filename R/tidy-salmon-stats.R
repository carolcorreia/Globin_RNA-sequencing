#####################
# Tidy Salmon stats #
#####################

# Author: Carolina N. Correia
# Date: August 3rd 2017

# Load packages
library(plyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(repurrrsive)

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")
getwd()

# Define variables for specific directories
summDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/summary_quant"
tablesDir <- c("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables")

# Create a vector with full paths to the salmon summary files
summ_files <- list.files(summDir, full.names = TRUE)

# Import and row-bind all files into a tibble
salmon_stats <- map_df(summ_files, ~ read_table2(.x))
View(salmon_stats)

# Rename library type
salmon_stats$Library_type %<>% 
    str_replace("IU", "Inward unstranded") %>% 
    str_replace("U", "Unstranded") %>% 
    str_replace("ISF", "Inward stranded forward")

# Add new column for plotting labels
salmon_stats$labels <- salmon_stats$File_name

# Correct plotting labels
salmon_stats$labels %<>%
    str_replace("_salmon_quant.log", "") %>% 
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



# Export renamed salmon stats
write_csv(salmon_stats,
          path = file.path(paste0(tablesDir, "/salmon_stats_tidy.csv")),
          col_names = TRUE)


# R session information
devtools::session_info()
