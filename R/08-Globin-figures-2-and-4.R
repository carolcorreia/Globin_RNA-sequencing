######################################
# Globin-derived transcripts RNA-seq #
#          Figures 2 and 4           #
######################################

# Author: Carolina N. Correia
# GitHub Repository DOI: 
# Date: September 9th 2017

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/R")

# Load previously plots
load("Plots-both-analyses.rda")

# Define figures directory
imgDir <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/Figures"

#############################
# 02 Load required packages #
#############################

# Load packages
library(devtools)
library(magrittr)
library(ggplot2)
library(ggjoy)
library(Cairo)
library(cowplot)

#install.packages(cowplot)

###############
# 03 Figure 2 #
###############

# Set grid
Fig2 <- plot_grid(jitter_plot,
                  jitter_plot2,
                  labels = c("A", "B"),
                  rel_heights = c(1, 0.8),
                  rel_widths = c(2, 0.8))

# Check plot
Fig2

# Export high quality PDF
ggsave("Figure_2-density.pdf",
       Fig2,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 12,
       units     = "in")

###############
# 04 Figure 4 #
###############

# Set grid
Fig4 <- plot_grid(joy_density,
                  joy_density2,
                  labels = c("A", "B"),
                  rel_heights = c(1, 0.8),
                  rel_widths = c(2, 1))

# Check plot
Fig4

# Export high quality PDF
ggsave("Figure_4-density.pdf",
       Fig4,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 12,
       units     = "in")

########################
# 05 Save .RData files #
########################

save.image(file = "Plots-both-analyses.rda")

#########################
# 06 Get R session info #
#########################

devtools::session_info()


