################################################################################
# ADNI LIPIDOMICS - VISUALIZATION OF AUTOSGI RESULTS
################################################################################
#
# Purpose: Generate detailed visualizations for selected clusters identified 
#          in AutoSGI analyses (ssGSEA, DPT, and pathway selection approaches)
#
# Study: Alzheimer's Disease Neuroimaging Initiative (ADNI)
#
# Visualizations generated:
#   - SGI tree plots showing hierarchical clustering and significance
#   - Overview plots with phenotype associations and feature heatmaps
#   - Feature abundance boxplots by cluster
#   - Phenotype distribution plots by cluster
#
# For each cluster of interest, creates three PDF files:
#   1. *_sgi_plots.pdf - Tree structure and overview
#   2. *_path_feat_plots.pdf - Feature abundances
#   3. *_pheno_plots.pdf - Phenotype distributions
#
################################################################################


################################################################################
# SETUP
## clear workspace and set directories ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
################################################################################

# Load required libraries
library(maplet)
library(rms)
library(stringr)
library(tidyverse)
library(magrittr)
library(openxlsx)
library(RColorBrewer)
library(sgi)
library(autosgi)

# Load custom functions
source("0_custom_functions.R")

## Define parameters ----

# Significance threshold for adjusted p-values
PADJ_THRESHOLD <- 0.05

# Output directory (from previous analysis run)
outdir <- paste0(Sys.Date(), "_adni-results/")

# Output file prefix
out_prefix <- "adni-lipidomics"


################################################################################
# DATA LOADING AND PREPARATION
################################################################################

## Load AutoSGI analysis results ----

# Load cleaned baseline data (from data cleaning step)
load(paste0(outdir, out_prefix, "-autosgi-data.rds"))
lipid_se_baseline <- se_bl; clinical_phenotypes <- adni_clinical
# Contains: lipid_se_baseline, colData_wide, clinical_phenotypes


# Load ssGSEA enrichment scores
load(paste0(outdir, out_prefix, "-autosgi-ssgsea-data.rds"))
# Contains: lipid_ssgsea_scores

# Load diffusion pseudotime scores
load(paste0(outdir, out_prefix, "-autosgi-dpt-data.rds"))
# Contains: pseudotime_matrix


## Prepare clinical phenotype data ----

# Select relevant clinical variables for plotting
# Exclude redundant diagnostic groups and MRI field variables
clinical_plot_data <- colData_wide %>% 
  select(
    -RID, -SC_DXGrp_bl, -SC_DXGrp_m12, -SC_DXGrp_m24, 
    -mri_bl, -mri_m12, -mri_m24
  )

# Set rownames for matching with lipid data
rownames(clinical_plot_data) <- colData_wide$RID


## Format clinical variable types ----

# Set ordinal class for diagnosis variables (severity ordering)
class(clinical_plot_data$DX_bl) <- "ordinal"
class(clinical_plot_data$DX_m12) <- "ordinal"
class(clinical_plot_data$DX_m24) <- "ordinal"

# Format MRI magnet strength as factor
clinical_plot_data$mag_bl <- as.factor(clinical_plot_data$mag_bl)
clinical_plot_data$mag_m12 <- as.factor(clinical_plot_data$mag_m12)
clinical_plot_data$mag_m24 <- as.factor(clinical_plot_data$mag_m24)

# Format demographic variables as factors
clinical_plot_data$PTGENDER <- as.factor(clinical_plot_data$PTGENDER)
clinical_plot_data$APOEGrp <- as.factor(clinical_plot_data$APOEGrp)


## Select phenotypes for visualization ----

# Focus on key clinical outcomes and covariates
# Includes diagnosis, hippocampal volume, tau biomarkers, and demographics
clinical_plot_data %<>% 
  select(
    # Diagnosis at each timepoint
    "DX_bl", "DX_m12", "DX_m24", 
    # Hippocampal volume at each timepoint
    "hip_vol_bl", "hip_vol_m12", "hip_vol_m24", 
    # Phospho-tau biomarkers at each timepoint
    "ptau_bl", "ptau_m12", "ptau_m24", 
    # Total tau biomarkers at each timepoint
    "tau_bl", "tau_m12", "tau_m24", 
    # Demographic and baseline variables
    "SC_Age", "PTGENDER", "PTEDUCAT", "APOEGrp", "bmi_bl"
  )


################################################################################
# VISUALIZE ssGSEA RESULTS
################################################################################
#
# Generates plots for enrichment-based clustering results
#
################################################################################

## Load ssGSEA cluster information ----

# Load cluster assignments from hierarchical selection
ssgsea_cluster_indices <- read.xlsx(
  paste0(outdir, out_prefix, "-ssgsea-hierarchical-selection-cluster-indices.xlsx")
)

# Path to AutoSGI association results
# Note: Use AutoSGI p-values (not SGI) because they're corrected across multiple trees
ssgsea_associations_file <- paste0(
  outdir, out_prefix, "-ssgsea-hierarchical-selection-association-results.xlsx"
)


## Define cluster of interest ----

# Cluster 18 identified as significant in ssGSEA analysis
ssgsea_cluster_id <- 18


## Generate ssGSEA visualizations ----

tree_specific_plots(
  method = "transform",                           # Using transformed enrichment scores
  cluster_of_interest = ssgsea_cluster_id, 
  main_df = lipid_ssgsea_scores,                 # ssGSEA enrichment matrix
  main_cl_idx = ssgsea_cluster_indices,          # Cluster membership assignments
  phenotype_data = clinical_plot_data,           # Clinical phenotypes
  padj_threshold = PADJ_THRESHOLD,               # Significance threshold
  asgi_file = ssgsea_associations_file,          # AutoSGI association results
  output_pdf_name = paste0(outdir, out_prefix, "-ssgsea-"), 
  pathway_column = "class",                       # Column with pathway annotations
  summarized_exp = lipid_se_baseline             # Original lipid data
)


################################################################################
# VISUALIZE DIFFUSION PSEUDOTIME (DPT) RESULTS
################################################################################
#
# Generates plots for pseudotime-based clustering results
#
################################################################################

## Load DPT cluster information ----

# Load cluster assignments from hierarchical selection
dpt_cluster_indices <- read.xlsx(
  paste0(outdir, out_prefix, "-dpt-hierarchical-selection-cluster-indices.xlsx")
)

# Path to AutoSGI association results
dpt_associations_file <- paste0(
  outdir, out_prefix, "-dpt-hierarchical-selection-association-results.xlsx"
)


## Define cluster of interest ----

# Cluster 17 identified as significant in DPT analysis
dpt_cluster_id <- 17


## Generate DPT visualizations ----

tree_specific_plots(
  method = "transform",                           # Using transformed pseudotime scores
  cluster_of_interest = dpt_cluster_id, 
  main_df = pseudotime_matrix,                   # DPT pseudotime matrix
  main_cl_idx = dpt_cluster_indices,             # Cluster membership assignments
  phenotype_data = clinical_plot_data,           # Clinical phenotypes
  padj_threshold = PADJ_THRESHOLD,               # Significance threshold
  asgi_file = dpt_associations_file,             # AutoSGI association results
  output_pdf_name = paste0(outdir, out_prefix, "-dpt-"), 
  pathway_column = "class",                       # Column with pathway annotations
  summarized_exp = lipid_se_baseline             # Original lipid data
)


################################################################################
# VISUALIZE PATHWAY SELECTION RESULTS
################################################################################
#
# Generates plots for direct pathway-level clustering results
#
################################################################################

## Load pathway selection information ----

# Path to pathway labels (mapping file)
pathway_labels_file <- paste0(
  outdir, out_prefix, "-pathway-selection-labels.xlsx"
)

# Path to AutoSGI association results
pathway_associations_file <- paste0(
  outdir, out_prefix, "-pathway-selection-association-results.xlsx"
)


## Define pathway of interest ----

# S1P (sphingosine-1-phosphate) pathway identified as significant
pathway_id <- "S1P"


## Generate pathway visualizations ----

tree_specific_plots(
  method = "pathway",                            # Using pathway-level selection
  cluster_of_interest = pathway_id,
  main_df = lipid_abundance_matrix,                  # Lipid abundance matrix
  main_cl_idx = NULL,            # Note: may need pathway-specific indices
  phenotype_data = clinical_plot_data,          # Clinical phenotypes
  padj_threshold = PADJ_THRESHOLD,              # Significance threshold
  asgi_file = pathway_associations_file,        # AutoSGI association results
  output_pdf_name = paste0(outdir, out_prefix, "-pathway-"),
  pathway_column = "class",                      # Column with pathway annotations
  summarized_exp = lipid_se_baseline,           # Original lipid data
  pathway_mapping = pathway_labels_file         # Pathway mapping file
)


################################################################################
# END OF SCRIPT
################################################################################
#
# Generated outputs for each cluster:
#   - {prefix}-{method}-{cluster}_sgi_plots.pdf
#   - {prefix}-{method}-{cluster}_path_feat_plots.pdf
#   - {prefix}-{method}-{cluster}_pheno_plots.pdf
#   - {prefix}-{method}-{cluster}_sgi_data.rds
#
################################################################################
