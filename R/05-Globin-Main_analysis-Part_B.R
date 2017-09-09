######################################
# Globin-derived transcripts RNA-seq #
#       Main Analysis - Part B       #
######################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: September 9th 2017

##################################
# 15 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Load previously saved data
load("Globin-Main_analysis.RData")

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
# 17 Keep genes with TPM >= 1 #
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
# 18 Summary statistics per treatment for each species #
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

########################################################
# 19 Plot: density of filtered gene counts per library #
########################################################

# Joyplot of density gene-level TPM after filtering
TPM_filt_all %>% 
    mutate(reverse_labels = fct_rev(labels)) %>% 
    ggplot(aes(x = log10(TPM + 1),
               y = reverse_labels)) +
        geom_joy(aes(fill = treatment), alpha = 0.5) +
        scale_fill_manual("Treatment",
                          values = c("#af8dc3", "#7fbf7b")) +
        theme_bw(base_size = 10) +
        ggtitle("Main analysis") +
        facet_wrap(~species, ncol = 2, scales = "free_y") +
        ylab("Density of gene-level TPM \nestimates per sample") +
        xlab(expression(paste(log[10], "(TPM + 1)"))) -> joy_density

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
# 20 Subset globin genes #
##########################

TPM_filt_all %>% 
    dplyr::filter(Gene_RefSeqID %in%
                      c("HBA", "HBA1", "HBA2", "HBB",
                        "LOC100737768", "LOC110259958",
                        "100036557", "100036558",
                        "100054109")) -> TPM_globins

# Check IDs
unique(TPM_globins$Gene_RefSeqID)

# Add gene symbols column for plotting
TPM_globins$gene_symbol <- TPM_globins$Gene_RefSeqID

TPM_globins$gene_symbol %<>% 
    str_replace("100036557","HBA1") %>% 
    str_replace("100036558","HBA2") %>% 
    str_replace("100054109","HBB") %>%
    str_replace("HBA$","HBA2") %>% 
    str_replace("LOC100737768","LOC100737768 (HBA)") %>%
    str_replace("LOC110259958","LOC110259958 (HBA)") %>%
    factor()

# Check factors
levels(TPM_globins$gene_symbol)
View(TPM_globins)

##################################################
# 21 Summary stats of globin genes per treatment #
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

# Pig Undepleted LOC100737768 (HBA)
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Porcine"
                    & gene_symbol == "LOC100737768 (HBA)")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "LOC100737768 (HBA)") %>% 
    dplyr::mutate(treatment = "Undepleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Pig Undepleted LOC110259958 (HBA)
TPM_globins %>%
    dplyr::filter(c(treatment == "Undepleted"
                    & species == "Porcine"
                    & gene_symbol == "LOC110259958 (HBA)")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "LOC110259958 (HBA)") %>% 
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

# Pig Globin depleted LOC100737768 (HBA)
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Porcine"
                    & gene_symbol == "LOC100737768 (HBA)")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "LOC100737768 (HBA)") %>% 
    dplyr::mutate(treatment = "Globin depleted") %>% 
    dplyr::bind_rows(globin_stats) -> globin_stats

# Pig Globin depleted LOC110259958 (HBA)
TPM_globins %>%
    dplyr::filter(c(treatment == "Globin depleted"
                    & species == "Porcine"
                    & gene_symbol == "LOC110259958 (HBA)")) %>%
    skim() %>%
    dplyr::select(c(var, stat, value)) %>% 
    dplyr::filter(c(var == "TPM" & stat != "hist")) %>% 
    dplyr::mutate(species = "Porcine") %>% 
    dplyr::mutate(gene_symbol = "LOC110259958 (HBA)") %>% 
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

#####################################
# 22 Correct globin summary factors #
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


# Treatment
globin_stats$treatment %<>% 
    factor() %>%
    forcats::fct_rev()

# Check species factors
levels(globin_stats$treatment)

###########################################
# 23 Tidy globins summary and export data #
###########################################

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
          file.path(paste0(tablesDir, "/globin-stats.csv")),
          col_names = TRUE)

##################################################
# 24 Plot: Distribution of globin gene-level TPM #
##################################################

jitter_plot <- ggplot(TPM_globins) +
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
    
    facet_grid(. ~ species, scales = "free_x") +
    theme_bw(base_size = 10) +
    ggtitle("Main analysis") +
    theme(axis.text.x = element_text(face = "italic",
                                     angle = 45,
                                     hjust = 1)) +
    xlab(NULL) +
    ylab(expression(paste(log[2], "(TPM + 1)"))) +
    labs(caption = expression(paste("Error bars represent the ",
                                    log[2], "(abs((mean + 1) ± SD))")))

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
# 25 Proportion of globin genes per sample #
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
    stringr::str_replace("Hsa_(S|P)\\d\\d_", "") %>% 
    stringr::str_replace("Ssc_\\d\\d_", "") %>%
    stringr::str_replace("U", "Undepleted") %>% 
    stringr::str_replace("D", "Globin depleted")

# Add species column
globin_proportion$species <- globin_proportion$labels

globin_proportion$species %<>% 
    stringr::str_replace("Bta_\\d\\d_(U|D)", "Bovine") %>% 
    stringr::str_replace("Eca_(T|S)\\d\\d_(U|D)", "Equine") %>% 
    stringr::str_replace("Hsa_(S|P)\\d\\d_(U|D)", "Human") %>% 
    stringr::str_replace("Ssc_\\d\\d_(U|D)", "Porcine")

# Export data
write_csv(globin_proportion,
          file.path(paste0(tablesDir, "/globin-proportion.csv")),
          col_names = TRUE)

################################################
# 26 Proportion of globin genes: summary stats #
################################################

globin_proportion %>% 
    dplyr::filter(species == "Human"
                  & treatment == "Undepleted") %>% 
    skim()

globin_proportion %>% 
    dplyr::filter(species == "Human"
                  & treatment == "Globin depleted") %>% 
    skim()

globin_proportion %>% 
    dplyr::filter(species == "Porcine"
                  & treatment == "Undepleted") %>% 
    skim()

globin_proportion %>% 
    dplyr::filter(species == "Porcine"
                  & treatment == "Globin depleted") %>% 
    skim()

globin_proportion %>% 
    dplyr::filter(species == "Equine"
                  & treatment == "Undepleted") %>% 
    skim()

globin_proportion %>% 
    dplyr::filter(species == "Bovine"
                  & treatment == "Undepleted") %>% 
    skim()

########################
# 27 Save .RData files #
########################

# Entire environment
save.image(file = "Globin-Main_analysis.RData")

# Plots only
save(jitter_plot, joy_density, file = "Plots-main-analysis.rda")

#########################
# 28 Get R session info #
#########################

devtools::session_info()

##################################
# Proceed to Additional analysis #
##################################

# File: 06-Globin-Main_analysis-Part_A.R



