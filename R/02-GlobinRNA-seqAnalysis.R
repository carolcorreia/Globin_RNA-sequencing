# HBA1, HBA2, and HBB transcripts RNA-seq Analysis - Part 2
# Author: Carolina C. Correia
# Date: August 17th 2017

##################################
# 15 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Load previously saved data
load("Globin-RNA-seqAnalysis.RData")

# Define variables for specific directories
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"
tablesDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/tables"

############################################
# 16 Load and/or install required packages #
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
# 17 Summary statistics per treatment for each species #
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
# 18 Plots: density of filtered gene counts per library #
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

##########################
# 19 Subset globin genes #
##########################

TPM_filt_all %>% 
    dplyr::filter(Gene_RefSeqID %in%
                      c("HBA","HBA1", "HBA2",
                        "HBB", "LOC100737768",
                        "100036557", "100036558",
                        "100054109")) -> TPM_globins

# Check IDs
unique(TPM_globins$Gene_RefSeqID)

# Add gene symbols column for plotting
TPM_globins$gene_symbol <- TPM_globins$Gene_RefSeqID

TPM_globins$gene_symbol %<>% 
    str_replace("LOC100737768","HBA1") %>%
    str_replace("100036557","HBA1") %>% 
    str_replace("100036558","HBA2") %>% 
    str_replace("100054109","HBB") %>%
    str_replace("HBA$","HBA2") %>% 
    factor()

# Check factors
levels(TPM_globins$gene_symbol)
View(TPM_globins)

##################################################
# 20 Summary stats of globin genes per treatment #
##################################################

# Cattle Undepleted HBB
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Bovine"
                    & gene_symbol == "HBB")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist"))  %>% 
    dplyr::mutate(species = "Bovine") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Undepleted") -> globin_stats

# Cattle Undepleted HBA2
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Bovine"
                    & gene_symbol == "HBA2")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist"))  %>% 
    dplyr::mutate(species = "Bovine") %>% 
    dplyr::mutate(gene_symbol = "HBA2") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Cattle Undepleted HBA1
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Bovine"
                    & gene_symbol == "HBA1")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Bovine") %>% 
    dplyr::mutate(gene_symbol = "HBA1") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Horse Undepleted HBB
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Equine"
                    & gene_symbol == "HBB")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Equine") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Horse Undepleted HBA2
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Equine"
                    & gene_symbol == "HBA2")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Equine") %>% 
    dplyr::mutate(gene_symbol = "HBA2") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Horse Undepleted HBA1
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Equine"
                    & gene_symbol == "HBA1")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Equine") %>% 
    dplyr::mutate(gene_symbol = "HBA1") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Pig Undepleted HBB
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Porcine"
                    & gene_symbol == "HBB")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Pig Undepleted HBA1
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Porcine"
                    & gene_symbol == "HBA1")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "HBA1") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Pig Globin depleted HBB
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Porcine"
                    & gene_symbol == "HBB")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Pig Globin depleted HBA1
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Porcine"
                    & gene_symbol == "HBA1")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "HBA1") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Human Undepleted HBB
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Human"
                    & gene_symbol == "HBB")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Human Undepleted HBA2
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Human"
                    & gene_symbol == "HBA2")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(gene_symbol = "HBA2") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Human Undepleted HBA1
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Human"
                    & gene_symbol == "HBA1")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(gene_symbol = "HBA1") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Human Globin depleted HBB
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Human"
                    & gene_symbol == "HBB")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(gene_symbol = "HBB") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Human Globin depleted HBA2
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Human"
                    & gene_symbol == "HBA2")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(gene_symbol = "HBA2") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Human Globin depleted HBA1
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Human"
                    & gene_symbol == "HBA1")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Human") %>% 
    dplyr::mutate(gene_symbol = "HBA1") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Check summary stats data frame
View(globin_stats)

#########################################################
# 21 Correct factors in globin summary stats data frame #
#########################################################

# Species
globin_stats$species %<>% 
    factor() %>%
    forcats::fct_inorder()

# Check species factors
levels(globin_stats$species)


# Gene symbol
globin_stats$gene_symbol %<>% 
    factor() %>%
    forcats::fct_inorder()

# Check gene symbol factors
levels(globin_stats$gene_symbol)


# Treatment
globin_stats$treatment %<>% 
    factor() %>%
    forcats::fct_rev()

# Check species factors
levels(globin_stats$treatment)

# Export

#################################################################
# 22 Spread mean and standard deviation of globins for plotting #
#################################################################

globin_stats %>% 
    dplyr::filter(stat == "mean") %>% 
    tidyr::spread(stat, value) -> mean_and_sd

globin_stats %>% 
    dplyr::filter(stat == "sd") %>% 
    tidyr::spread(stat, value) %>% 
    dplyr::right_join(mean_and_sd) -> mean_and_sd

# Check data frame
mean_and_sd

mean_and_sd %>% 
dplyr::filter(treatment == "Undepleted") -> mean_and_sd_U
#########################################################
#  Jitter plot #
#########################################################



jitter_plot <- ggplot(TPM_globins) +
    geom_jitter(aes(gene_symbol, log2(TPM + 1),
                    colour = treatment),
                size = 3,
                alpha = 0.7) +
    scale_colour_manual("Treatment",
                        values = c("#af8dc3", "#7fbf7b")) +
    geom_errorbar(data = mean_and_sd,
                  aes(x = gene_symbol, ymin = log2(abs((mean + 1) - sd)),
                      ymax= log2(abs((mean + 1) + sd))),
                  colour = "#414545",
                  width = 0.3) +
    #scale_errorbar_manual(expression(paste(log[2], "(test)"))) +
    geom_point(data = mean_and_sd, aes(gene_symbol, log2(mean + 1),
                  shape = treatment),
               colour = "black",
              size = 3) +
    scale_shape_manual(expression(paste(log[2], "(mean + 1)")),
                       values = c(17, 15)) +
    
    facet_grid(. ~ species, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(face = "italic")) +
    xlab(NULL) +
    ylab(expression(paste(log[2], "(TPM + 1)"))) +
    labs(title = "Globin genes", 
     subtitle = "test", 
     caption = expression(paste("Error bars represent the ",
                                log[2], "(abs((mean + 1) ± SD))")))

# Change title and subtitle
jitter_plot


ggsave(".svg",
       plot      = plot,
       limitsize = FALSE,
       dpi       = 300,
       path      = img_dir)


#######################
#  Save .RData file #
#######################

save.image(file = "Globin-RNA-seqAnalysis.RData")

#########################
#  Get R session info #
#########################

devtools::session_info()



