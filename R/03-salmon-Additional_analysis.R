#####################################################
#    Globin-derived transcripts RNA-seq Analysis    #
# Salmon mapping statistics for additional analysis #
#####################################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: August 29th 2017

# Load packages
library(plyr)
library(tidyverse)
library(magrittr)
library(stringr)
library(skimr)
library(devtools)

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/extra_horse_summary")

# Define variables for specific directories
summDir <- getwd()
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables"

# Create a vector with full paths to the salmon summary files
summ_file <- list.files(summDir, full.names = TRUE)

# Import and row-bind all files into a tibble
salmon_stats <- read.table(summ_file, header = TRUE)
View(salmon_stats)

# Rename columns
salmon_stats %<>% 
    dplyr::rename(Observed_fragments = Total_fragments) %>% 
    dplyr::rename(Mapped_fragments = Total_reads) %>% 
    dplyr::rename(`Mapping_rate(%)` = `Mapping_rate...`) %>% 
    as.tibble()

# Remove % symbol and convert mapping rate to numeric
salmon_stats$`Mapping_rate(%)` %<>%
    str_replace("%", "") %>% 
    as.numeric()

# Rename library type
salmon_stats$Library_type %<>% 
    str_replace("U", "Unstranded")

# Add new column for plotting labels
salmon_stats$labels <- salmon_stats$File_name

# Correct plotting labels
salmon_stats$labels %<>%
    str_replace("_salmon_quant.log", "") %>% 
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
    stringr::str_replace("Eca_(T|S)\\d\\d_", "") %>% 
    stringr::str_replace("U", "Undepleted")

# Add species column
salmon_stats$species <- salmon_stats$labels

salmon_stats$species %<>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_(U|D)", "Equine")

# Check data frame
View(salmon_stats)

# Export complete stats
write_csv(salmon_stats,
          path = file.path(paste0(tablesDir, "/horse-extra-salmon-stats.csv")),
          col_names = TRUE)

# Calculate summary stats
vars_to_keep <- c("Mapping_rate(%)", "Observed_fragments", "Mapped_fragments")

salmon_stats %>%
    skim() %>%
    dplyr::filter(var %in% vars_to_keep) %>% 
    dplyr::filter(stat == "mean") %>% 
    dplyr::mutate(stat = paste(stat, level, sep = '_')) %>% 
    dplyr::select(-c(type, level)) %>% 
    dplyr::mutate(var = paste(var, stat, sep = '_')) %>% 
    dplyr::select(-stat) %>% 
    tidyr::spread(var, value) %>% 
    dplyr::mutate(species = "Equine") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::mutate(library_type = "Unstranded") -> salmon_summary

# Reorder columns in summary data frame
salmon_summary %<>% 
    dplyr::rename(Observed_fragments_mean = Observed_fragments_mean_.all,
                  Mapped_fragments_mean = Mapped_fragments_mean_.all,
                  `Mapping_rate(%)_mean` = `Mapping_rate(%)_mean_.all`) %>%
    dplyr::select(species, treatment,
                  library_type, Observed_fragments_mean,
                  everything())

# Check data frame
View(salmon_summary)

# Export summary stats
write_csv(salmon_summary,
          path = file.path(paste0(tablesDir, "/extra-horse-mean.csv")),
          col_names = TRUE)

# R session information
devtools::session_info()
