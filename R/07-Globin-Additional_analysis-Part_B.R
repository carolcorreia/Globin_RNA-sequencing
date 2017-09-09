######################################
# Globin-derived transcripts RNA-seq #
#   Additional Analysis - Part B     #
######################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: September 9th 2017

##################################
# 13 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Load previously saved data
load("Globin-Additional_analysis.RData")

# Define variables for specific directories
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables"

############################################
# 14 Load and/or install required packages #
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

#######################################
# 15 Summary statistics per treatment #
#######################################

# Calculate summary stats
horse_filt %>%
    skim() -> horse_stats


# Visualise stats
skim_print(horse_stats)

# Export stats
write_csv(horse_stats,
          file.path(paste0(tablesDir, "/", "horse_extra_stats.csv")),
          col_names = TRUE)

#########################################################
# 16 Plot: density of filtered gene counts per library #
#########################################################

# Joyplot of density gene-level TPM after filtering
horse_filt %>% 
    mutate(reverse_labels = fct_rev(labels)) %>% 
    ggplot(aes(x = log10(TPM + 1),
               y = reverse_labels)) +
        geom_joy(aes(fill = treatment), alpha = 0.5) +
        scale_fill_manual("Treatment",
                          values = c("#af8dc3", "#7fbf7b")) +
        theme_bw(base_size = 10) +
        ggtitle("Additional analysis") +
        ylab("Density of gene-level TPM \nestimates per sample") +
        xlab(expression(paste(log[10], "(TPM + 1)"))) -> joy_density2

# Check plot
joy_density2

# Export high quality PDF
ggsave("horse-extra-joy-density.pdf",
       joy_density2,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 4,
       units     = "in")

######################
# 17 Subset HBB gene #
######################

horse_filt %>% 
    dplyr::filter(Gene_RefSeqID %in%
                      c("HBA", "HBA1", "HBA2", "HBB")) -> TPM_globins

# Check IDs
unique(TPM_globins$Gene_RefSeqID)

# Add gene symbols column for plotting
TPM_globins$gene_symbol <- TPM_globins$Gene_RefSeqID

TPM_globins$gene_symbol %<>% 
    factor()

# Check data frame
View(TPM_globins)

###########################
# 18 Summary stats of HBB #
###########################

TPM_globins %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Equine") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    as.tibble() -> globin_stats

# Check summary stats data frame
View(globin_stats)

#######################################
# 19 Tidy HBB summary and export data #
#######################################

globin_stats %>% 
    dplyr::filter(stat == "mean") %>% 
    tidyr::spread(stat, value) -> globin_tidy

globin_stats %>% 
    dplyr::filter(stat == "median") %>% 
    tidyr::spread(stat, value) %>% 
    dplyr::right_join(globin_tidy) -> globin_tidy

globin_stats %>% 
    dplyr::filter(stat == "sd") %>% 
    tidyr::spread(stat, value) %>% 
    dplyr::right_join(globin_tidy) -> globin_tidy

globin_stats %>% 
    dplyr::filter(stat == "n") %>% 
    tidyr::spread(stat, value) %>% 
    dplyr::right_join(globin_tidy) -> globin_tidy

globin_stats %>% 
    dplyr::filter(stat == "min") %>% 
    tidyr::spread(stat, value) %>% 
    dplyr::right_join(globin_tidy) -> globin_tidy

globin_stats %>% 
    dplyr::filter(stat == "max") %>% 
    tidyr::spread(stat, value) %>% 
    dplyr::right_join(globin_tidy) -> globin_tidy

# Check data frame
globin_tidy

# Export tidy globin summary stats
write_csv(globin_tidy,
          file.path(paste0(tablesDir, "/horse-extra-globin-stats.csv")),
          col_names = TRUE)

###############################################
# 20 Plot: Distribution of HBB gene-level TPM #
###############################################

jitter_plot2 <- ggplot(TPM_globins) +
    geom_jitter(aes(gene_symbol, log2(TPM + 1),
                    colour = treatment),
                size = 3,
                alpha = 0.7) +
    scale_colour_manual("Treatment",
                        values = c("#af8dc3", "#7fbf7b")) +
    geom_errorbar(data = globin_tidy,
                  aes(x = gene_symbol, ymin = log2(abs((mean + 1) - sd)),
                      ymax= log2(abs((mean + 1) + sd))),
                  colour = "#414545",
                  width = 0.3) +
    geom_point(data = globin_tidy, aes(gene_symbol, log2(mean + 1),
                  shape = treatment),
               colour = "black",
              size = 3) +
    scale_shape_manual(expression(paste(log[2], "(mean + 1)")),
                       values = c(17, 15)) +
    theme_bw(base_size = 10) +
    ggtitle("Additional analysis") +
    theme(axis.text.x = element_text(face = "italic",
                                     angle = 45,
                                     hjust = 1)) +
    xlab(NULL) +
    ylab(expression(paste(log[2], "(TPM + 1)"))) +
    labs(caption = expression(paste("Error bars represent the ",
                                    log[2], "(abs((mean + 1) ± SD))")))

# Check plot
jitter_plot2

# Export high quality PDF
ggsave("horse-extra-jitter.pdf",
       jitter_plot2,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 4,
       width     = 4,
       units     = "in")

###################################
# 21 Proportion of HBB per sample #
###################################

# Get total TPM per sample
horse_filt %>% 
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

# Export data
write_csv(globin_proportion,
          file.path(paste0(tablesDir, "/horse-extra-HBB-prop.csv")),
          col_names = TRUE)

#######################################
# 22 Proportion of HBB: summary stats #
#######################################

globin_proportion %>% 
    skim() -> globin_summary

# Chech data frame
skim_print(globin_summary)

########################
# 23 Save .RData files #
########################

# Entire environment
save.image(file = "Globin-Additional_analysis.RData")

# Only plots
save(joy_density2, jitter_plot2, file = "Plots-both-analyses.rda")

#########################
# 24 Get R session info #
#########################

devtools::session_info()



