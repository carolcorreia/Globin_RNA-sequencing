# HBA1, HBA2, and HBB transcripts RNA-seq Analysis - Part 2
# Author: Carolina C. Correia
# Date: August 15th 2017

##################################
# 14 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")
getwd()

# Load previously saved data
load("Globin-RNA-seqAnalysis.RData")

# Define variables for specific directories
workDir <- getwd()
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables"

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
library(ggforce)
library(waffle)
library(skimr)

# Uncomment functions below to install packages in case you don't have them

#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("ggjoy")
#install.packages("ggforce")
#library(devtools)
#devtools::install_github("hrbrmstr/waffle")
#devtools::install_github("hadley/colformat")
#devtools::install_github("ropenscilabs/skimr")

########################################################
# 16 Summary statistics per treatment for each species #
########################################################

# Calculate summary stats
TPM_filt_all %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Human")) %>%
    skim() -> human_U

TPM_filt_all %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Human")) %>%
    skim() -> human_D

TPM_filt_all %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Porcine")) %>%
    skim() -> pig_U

TPM_filt_all %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Porcine")) %>%
    skim() -> pig_D

TPM_filt_all %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Equine")) %>%
    skim() -> horse_U

TPM_filt_all %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Bovine")) %>%
    skim() -> cattle_U

# Visualise stats
skim_print(human_U)
skim_print(human_D)
skim_print(pig_U)
skim_print(pig_D)
skim_print(horse_U)
skim_print(cattle_U)

# Export stats for all subsets
vars_summary <- list(human_U, human_D, pig_U, pig_D, horse_U, cattle_U)
files_summ <- paste0(c("human_U", "human_D",
                       "pig_U", "pig_D",
                       "horse_U", "cattle_U"), "_stats.csv")

path_summ <- file.path(paste0(tablesDir, "/", files_summ))

purrr::pwalk(list(vars_summary, path_summ),
             write_csv,
             col_names = TRUE)

#########################################################
# 17 Plots: density of filtered gene counts per library #
#########################################################

# Joyplot of density gene-level TPM after filtering
TPM_filt_all %>% 
    mutate(reverse_labels = fct_rev(labels)) %>% 
    ggplot(aes(x = log10(TPM + 1),
               y = reverse_labels)) +
        geom_joy(aes(fill = treatment), alpha = 0.5) +
        scale_fill_manual("Treatment",
                          values = c("#af8dc3", "#7fbf7b")) +
        theme_bw() +
        facet_wrap(~species, ncol = 2, scales = "free_y") +
        ylab("Density of gene-level TPM \nestimates per sample") +
        xlab(expression(paste(log[10], "(TPM + 1)"))) -> joy_density

# Standard density plot of gene-level TPM after filtering
density_plot <- ggplot(TPM_filt_all) +
                geom_density(aes(log10(TPM + 1),
                                 group = labels,
                                 colour = treatment),
                             alpha = 0.5, show.legend = FALSE) +
                stat_density(aes(x = log10(TPM + 1),
                                 colour = treatment),
                             geom = "line", position = "identity") +
                scale_colour_manual("Treatment",
                                    values = c("#af8dc3", "#7fbf7b")) +
                guides(colour = guide_legend(override.aes = list(size = 3))) +
                facet_grid(. ~ species) +
                theme_bw() +
                ylab("Density of gene-level TPM \nestimates per sample") +
                xlab(expression(paste(log[10], "(TPM + 1)")))

# Export high quality image for both plots
vars_den <- list(joy_density, density_plot)

# PNG
png_den <- paste0(c("joy_density", "density_plot"), ".png")
purrr::pwalk(list(png_den, vars_den),
             ggsave,
             path      = imgDir,
             device    = "png",
             height    = 30,
             width     = 30,
             units     = "cm",
             limitsize = FALSE,
             dpi       = 600)

# SVG
svg_den <- paste0(c("joy_density", "density_plot"), ".svg")
purrr::pwalk(list(svg_den, vars_den),
             ggsave,
             path      = imgDir,
             device    = "svg",
             limitsize = FALSE,
             dpi       = 600)




#######################
#  Save .RData file #
#######################

save.image(file = "Globin-RNA-seqAnalysis.RData")

#########################
#  Get R session info #
#########################

devtools::session_info()



