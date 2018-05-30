######################################
# Globin-derived transcripts RNA-seq #
#       2018 Analysis - Part B       #
######################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: May 29th 2018

##################################
# 14 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Load previously saved data
load("Globin-2018_analysis.RData")

# Define variables for specific directories
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures_2018"
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables_2018"

############################################
# 15 Load and/or install required packages #
############################################

# Load packages
library(plyr)
library(tidyverse)
library(devtools)
library(magrittr)
library(stringr)
library(forcats)
library(ggjoy)
library(Cairo)
library(skimr)

# Uncomment functions below to install packages in case you don't have them

#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("ggjoy")
#install.packages("Cairo")
#library(devtools)
#devtools::install_github("hadley/colformat")
#devtools::install_github("ropenscilabs/skimr")

###############################
# 16 Keep genes with TPM >= 1 #
###############################

TPM_filt_all %>% 
    dplyr::group_by(labels) %>% 
    dplyr::count(TPM >= 1) %>% 
    dplyr::filter(`TPM >= 1` == TRUE) %>% 
    write_csv(file.path(paste0(tablesDir, "/Expressed-genes-sample.csv")),
              col_names = TRUE)

TPM_filt_all %<>% 
    dplyr::group_by(labels) %>% 
    dplyr::filter(TPM >= 1)

########################################################
# 17 Summary statistics per treatment for each species #
########################################################

# Calculate summary stats
TPM_filt_all %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Equine")) %>%
    skim() -> horse_U

TPM_filt_all %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Bovine")) %>%
    skim() -> cattle_U

# Visualise stats
skim(horse_U)
skim(cattle_U)

# Export stats for all subsets
vars_summary <- list(horse_U, cattle_U)
files_summ <- paste0(c("horse_U", "cattle_U"), "_stats.csv")

path_summ <- file.path(paste0(tablesDir, "/", files_summ))

purrr::pwalk(list(vars_summary, path_summ),
             write_csv,
             col_names = TRUE)

########################################################
# 18 Plot: density of filtered gene counts per library #
########################################################

# Joyplot of density gene-level TPM after filtering
TPM_filt_all %>% 
    mutate(reverse_labels = fct_rev(labels)) %>% 
    ggplot(aes(x = log10(TPM),
               y = reverse_labels)) +
        geom_joy(aes(fill = treatment), alpha = 0.5) +
        scale_fill_manual("Treatment",
                          values = c("#af8dc3", "#7fbf7b")) +
        theme_bw(base_size = 10) +
        ggtitle("2018 analysis") +
        facet_wrap(~species, ncol = 2, scales = "free_y") +
        ylab("Density of gene-level TPM \nestimates per sample") +
        xlab(expression(paste(log[10], "TPM"))) -> joy_density

# Check plot
joy_density

# Export high quality PDF
ggsave("joy_density.pdf",
       joy_density,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

##########################
# 19 Subset globin genes #
##########################

TPM_filt_all %>% 
    dplyr::filter(Gene_RefSeqID %in%
                      c("HBA", "HBA1", "HBA2", "HBB")) -> TPM_globins

# Check IDs
unique(TPM_globins$Gene_RefSeqID)

# Add gene symbols column for plotting
TPM_globins$gene_symbol <- TPM_globins$Gene_RefSeqID

TPM_globins %<>% 
    mutate(gene_symbol = ifelse(species == "Equine" &
                                    gene_symbol == "HBA",
                                "HBA (HBA1)", gene_symbol)) %>% 
    mutate(gene_symbol = ifelse(species == "Bovine" &
                                    gene_symbol == "HBA",
                                "HBA (HBA2)", gene_symbol))

# Create factors
TPM_globins$gene_symbol %<>% 
    factor(levels = c("HBA (HBA1)", "HBA2",
           "HBA1", "HBA (HBA2)", "HBB"))

# Check factors
levels(TPM_globins$gene_symbol)
View(TPM_globins)

##################################################
# 20 Summary stats of globin genes per treatment #
##################################################

# Cattle
TPM_globins %>%
    dplyr::filter(species == "Bovine") %>% 
    group_by(gene_symbol) %>%
    skim_to_wide() %>% 
    dplyr::mutate(species = "Bovine") -> globin_stats 

# Horse
TPM_globins %>%
    dplyr::filter(species == "Equine") %>% 
    group_by(gene_symbol) %>% 
    skim_to_wide() %>% 
    dplyr::mutate(species = "Equine") %>%
    dplyr::bind_rows(globin_stats) -> globin_stats

# Check summary stats data frame
View(globin_stats)

#####################################
# 21 Correct globin summary factors #
#####################################

# Species
globin_stats$species %<>% 
    factor() %>%
    forcats::fct_inorder()

# Check species factors
levels(globin_stats$species)


# Gene symbol
globin_stats$gene_symbol %<>% 
    factor() %>%
    forcats::fct_inorder() %>% 
    forcats::fct_relevel("HBB", after = Inf)

# Check gene symbol factors
levels(globin_stats$gene_symbol)

###########################################
# 22 Tidy globins summary and export data #
###########################################

globin_stats %>% 
    dplyr::filter(variable == "TPM") %>% 
    dplyr::select(-c(type, hist)) %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::mutate(mean = as.numeric(mean),
                  sd = as.numeric(sd)) -> globin_tidy

# Check data frame
globin_tidy

# Export tidy globin summary stats
write_csv(globin_tidy,
          file.path(paste0(tablesDir, "/globin-stats.csv")),
          col_names = TRUE)

##################################################
# 23 Plot: Distribution of globin gene-level TPM #
##################################################

jitter_plot <- ggplot(TPM_globins) +
    geom_jitter(aes(gene_symbol, log2(TPM),
                    colour = treatment),
                size = 3,
                alpha = 0.7) +
    scale_colour_manual("Treatment",
                        values = c("#af8dc3", "#7fbf7b")) +
    geom_errorbar(data = globin_tidy,
                  aes(x = gene_symbol, ymin = log2(abs(mean - sd)),
                      ymax= log2(abs(mean + sd))),
                  colour = "#414545",
                  width = 0.3) +
    geom_point(data = globin_tidy, aes(gene_symbol, log2(mean),
                  shape = treatment),
               colour = "black",
              size = 3) +
    scale_shape_manual(expression(paste(log[2], "mean")),
                       values = c(17, 15)) +
    
    facet_grid(. ~ species, scales = "free_x") +
    theme_bw(base_size = 10) +
    ggtitle("2018 analysis") +
    theme(axis.text.x = element_text(face = "italic",
                                     angle = 45,
                                     hjust = 1)) +
    xlab(NULL) +
    ylab(expression(paste(log[2], "TPM"))) +
    labs(caption = expression(paste("Error bars represent the ",
                                    log[2], "(abs(mean ± SD))")))

# Check plot
jitter_plot

# Export high quality PDF
ggsave("jitter_plot.pdf",
       jitter_plot,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 7,
       width     = 10,
       units     = "in")

############################################
# 24 Proportion of globin genes per sample #
############################################

# Get total TPM per sample
TPM_filt_all %>% 
    dplyr::group_by(labels) %>% 
    dplyr::summarise(Total_TPM_all_genes = sum(TPM)) -> globin_proportion

# Get total globins TPMs per sample
TPM_globins %>% 
    dplyr::group_by(labels) %>% 
    dplyr::summarise(Total_TPM_globins = sum(TPM)) %>% 
    dplyr::right_join(globin_proportion) -> globin_proportion

# Order columns and calculate proportion (%)
globin_proportion %<>% 
    dplyr::select(labels, Total_TPM_all_genes, Total_TPM_globins) %>% 
    dplyr::mutate(Percent_globins =
                      Total_TPM_globins / Total_TPM_all_genes * 100)

# Check data frame
View(globin_proportion)

# Add treatment column and move it
globin_proportion$treatment <- globin_proportion$labels

globin_proportion %<>% 
    dplyr::select(labels, treatment, everything())

globin_proportion$treatment %<>% 
    stringr::str_replace("Bta_\\d\\d_", "") %>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_", "") %>% 
    stringr::str_replace("U", "Undepleted")

# Add species column
globin_proportion$species <- globin_proportion$labels

globin_proportion$species %<>% 
    stringr::str_replace("Bta_\\d\\d_(U|D)", "Bovine") %>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_(U|D)", "Equine")

# Export data
write_csv(globin_proportion,
          file.path(paste0(tablesDir, "/globin-proportion.csv")),
          col_names = TRUE)

################################################
# 25 Proportion of globin genes: summary stats #
################################################

globin_proportion %>% 
    dplyr::filter(species == "Equine"
                  & treatment == "Undepleted") %>% 
    skim()

globin_proportion %>% 
    dplyr::filter(species == "Bovine"
                  & treatment == "Undepleted") %>% 
    skim()

########################
# 26 Save .RData files #
########################

# Entire environment
save.image(file = "Globin-2018_analysis.RData")

# Selected variables
save(jitter_plot,
     joy_density,
     TPM_filt_all,
     file = "Plots-2018-analysis.rda")

#########################
# 27 Get R session info #
#########################

devtools::session_info()





