################################################################################
# ROSMAP METABOLOMICS - VISUALIZATION OF AUTOSGI RESULTS
################################################################################
#
# Purpose: Generate detailed visualizations for selected clusters identified 
#          in AutoSGI hierarchical selection analysis of ROS/MAP brain metabolomics
#
# Study: Religious Orders Study and Memory and Aging Project (ROS/MAP)
#
# Last modified: Richa Batra
# Date: 2026-02-07
#
# Visualizations generated:
#   - SGI tree plots showing hierarchical clustering and significance
#   - Overview plots with phenotype associations and feature heatmaps
#   - Feature abundance boxplots by cluster
#   - Phenotype distribution plots by cluster
#
# For each cluster of interest, creates three PDF files:
#   1. *_sgi_plots.pdf - Tree structure and overview
#   2. *_path_feat_plots.pdf - Metabolite abundances
#   3. *_tree_plots.pdf - Phenotype distributions
#
################################################################################


################################################################################
# SETUP
################################################################################

## Clear workspace and set working directory ----

# Remove all objects from workspace for clean analysis
rm(list = ls())

# Set working directory to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


## Load required libraries ----

library(maplet)
library(sgi)
library(autosgi)
library(rms)
library(stringr)
library(tidyverse)
library(magrittr)
library(openxlsx)
library(RColorBrewer)

# Load custom functions
source("0_custom_functions.R")


## Define parameters ----

# Significance threshold for adjusted p-values
PADJ_THRESHOLD <- 0.05

# Output directory (from previous analysis run)
outdir <- paste0(Sys.Date(), "_rosmap-results/")

# Output file prefix
out_prefix <- "rosmap-metabo"


################################################################################
# DATA LOADING AND PREPARATION
################################################################################

## Load AutoSGI analysis results ----

# Load cleaned data from hierarchical selection analysis
load(paste0(outdir, out_prefix, "-autosgi-data.rds"))
# Contains: rosmap_data, metabolite_matrix, phenotype_data


## Format clinical variable types ----

# Set ordinal class for variables with natural severity ordering
class(phenotype_data$cogdx) <- "ordinal"    # Cognitive diagnosis severity
class(phenotype_data$braaksc) <- "ordinal"  # Braak staging (tau pathology, 0-VI)
class(phenotype_data$ceradsc) <- "ordinal"  # CERAD score (amyloid pathology severity)


## Load cluster information ----

# Load cluster assignments from hierarchical selection
cluster_indices <- read.xlsx(
  paste0(outdir, out_prefix, "-hierarchical-selection-cluster-indices.xlsx")
)

# Path to AutoSGI association results
# Note: Use AutoSGI p-values (not SGI) because they're corrected across multiple trees
associations_file <- paste0(
  outdir, out_prefix, "-hierarchical-selection-association-results.xlsx"
)


################################################################################
# VISUALIZE HIERARCHICAL SELECTION RESULTS
################################################################################
#
# Generates plots for metabolite clusters identified through hierarchical
# clustering and supervised group identification
#
################################################################################

## Define cluster of interest ----

# Cluster 358 identified as significantly associated with clinical outcomes
# in hierarchical selection analysis
cluster_id <- 358


## Generate cluster-specific visualizations ----

tree_specific_plots(
  method = "hc_sel",                              # Using hierarchical clustering selection
  cluster_of_interest = cluster_id, 
  main_df = metabolite_matrix,                    # Metabolite abundance matrix
  main_cl_idx = cluster_indices,                  # Cluster membership assignments
  phenotype_data = phenotype_data,                # Clinical phenotypes
  padj_threshold = PADJ_THRESHOLD,                # Significance threshold
  asgi_file = associations_file,                  # AutoSGI association results
  output_pdf_name = paste0(outdir, out_prefix, "-hierarchical-selection-"), 
  pathway_column=NULL,
  summarized_exp = rosmap_data,                    # Original SummarizedExperiment
  pathway_mapping = NULL
)


################################################################################
# END OF SCRIPT
################################################################################
#
# Generated outputs for cluster 358:
#   - rosmap-metabo-hierarchical-selection-358_sgi_plots.pdf
#   - rosmap-metabo-hierarchical-selection-358_path_feat_plots.pdf
#   - rosmap-metabo-hierarchical-selection-358_tree_plots.pdf
#   - rosmap-metabo-hierarchical-selection-358_sgi_data.rds
#
################################################################################
