###############################################
# Globin-derived transcripts RNA-seq Analysis #
# Salmon mapping statistics for main analysis #
###############################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: June 20th 2018

# Load packages
library(plyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(skimr)
library(devtools)

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/summary_quant")

# Define variables for specific directories
summDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/summary_quant/2018"
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables_2018"

# Create a vector with full paths to the salmon summary files
summ_files <- c(list.files(summDir, full.names = TRUE))

# Import and row-bind all files into a tibble
salmon_stats <- map_df(summ_files, ~ read_table2(.x))
View(salmon_stats)

# Rename columns
salmon_stats %<>% 
    dplyr::rename(Observed_fragments = Total_fragments) %>% 
    dplyr::rename(Mapped_fragments = Total_reads)

# Remove % symbol and convert mapping rate to numeric
salmon_stats$`Mapping_rate(%)` %<>%
    str_replace("%", "") %>% 
    as.numeric()

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


# Check data frame
View(salmon_stats)

# Convert labels into factors and order them
salmon_stats$labels %<>%
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
levels(salmon_stats$labels)

# Move labels column to the beginning and order rows by labels
salmon_stats %<>% 
    dplyr::select(labels, everything()) %>% 
    dplyr::arrange(labels)

# Add treatment column and move it
salmon_stats$treatment <- salmon_stats$labels

salmon_stats %<>% 
    dplyr::select(labels, treatment, everything())

salmon_stats$treatment %<>% 
    stringr::str_replace("Bta_\\d\\d_", "") %>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_", "") %>% 
    stringr::str_replace("Hsa_(S|P)\\d\\d_", "") %>% 
    stringr::str_replace("Ssc_\\d\\d_", "") %>%
    stringr::str_replace("U", "Undepleted") %>% 
    stringr::str_replace("D", "Globin depleted")

# Add species column
salmon_stats$species <- salmon_stats$labels

salmon_stats$species %<>% 
    stringr::str_replace("Bta_\\d\\d_(U|D)", "Bovine") %>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_(U|D)", "Equine") %>% 
    stringr::str_replace("Hsa_(S|P)\\d\\d_(U|D)", "Human") %>% 
    stringr::str_replace("Ssc_\\d\\d_(U|D)", "Porcine")

# Check data frame
View(salmon_stats)

# Export complete stats
write_csv(salmon_stats,
          path = file.path(paste0(tablesDir, "/salmon_stats_tidy.csv")),
          col_names = TRUE)

# Calculate summary stats
vars_to_keep <- c("Mapping_rate(%)", "Observed_fragments", "Mapped_fragments")

salmon_stats %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Bovine")) %>%
    skim() %>%
    dplyr::filter(variable %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::select(-c(type, level, stat, formatted)) %>% 
    tidyr::spread(variable, value) %>% 
    dplyr::mutate(species = "Bovine") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::mutate(library_type = "Inward stranded forward") -> salmon_summary
    
salmon_stats %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Equine")) %>%
    skim() %>%
    dplyr::filter(variable %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::select(-c(type, level, stat, formatted)) %>% 
    tidyr::spread(variable, value) %>% 
    dplyr::mutate(species = "Equine") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::mutate(library_type = "Unstranded") %>% 
    dplyr::bind_rows(salmon_summary) -> salmon_summary

salmon_stats %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Porcine")) %>%
    skim() %>%
    dplyr::filter(variable %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::select(-c(type, level, stat, formatted)) %>% 
    tidyr::spread(variable, value) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::mutate(library_type = "Inward unstranded") %>% 
    dplyr::bind_rows(salmon_summary) -> salmon_summary

salmon_stats %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Porcine")) %>%
    skim() %>%
    dplyr::filter(variable %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::select(-c(type, level, stat, formatted)) %>% 
    tidyr::spread(variable, value) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::mutate(library_type = "Inward unstranded") %>% 
    dplyr::bind_rows(salmon_summary) -> salmon_summary

salmon_stats %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Human")) %>%
    skim() %>%
    dplyr::filter(variable %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::select(-c(type, level, stat, formatted)) %>% 
    tidyr::spread(variable, value) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::mutate(library_type = "Inward unstranded") %>% 
    dplyr::bind_rows(salmon_summary) -> salmon_summary

salmon_stats %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Human")) %>%
    skim() %>%
    dplyr::filter(variable %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::select(-c(type, level, stat, formatted)) %>% 
    tidyr::spread(variable, value) %>%  
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::mutate(library_type = "Inward unstranded") %>% 
    dplyr::bind_rows(salmon_summary) -> salmon_summary

# Reorder columns in summary data frame
salmon_summary %<>% 
    dplyr::rename(Observed_fragments_mean = Observed_fragments,
                  Mapped_fragments_mean = Mapped_fragments,
                  `Mapping_rate(%)_mean` = `Mapping_rate(%)`) %>%
    dplyr::select(species, treatment,
                  library_type, Observed_fragments_mean,
                  everything())
    
# Check data frame
View(salmon_summary)

# Export summary stats
write_csv(salmon_summary,
          path = file.path(paste0(tablesDir, "/salmon_mean.csv")),
          col_names = TRUE)

# R session information
devtools::session_info()
