################################################################################
# ADNI LIPIDOMICS - AUTOSGI CASE STUDY (not reported)
################################################################################
#
# Purpose: Analyze ADNI lipidomics data using three complementary approaches:
#          1. ssGSEA - Pathway enrichment-based clustering
#          2. DPT - Diffusion pseudotime-based clustering  
#          3. Pathway Selection - Direct pathway-level clustering
#
# Study: Alzheimer's Disease Neuroimaging Initiative (ADNI)
# Data: Baker Lipidomics Lab longitudinal measurements from LONI
# Input file is a summarized experiment in an Excel file. This data corresponds to the 
# data from the publication PMID: 40592256. Data is available through LONI and sage/synapse. 
# More information on data availablility in the manuscript. 
#
#
# Analysis Rationale:
#   - ssGSEA: Identifies lipid class-level enrichment patterns
#   - DPT: Captures temporal/trajectory patterns in lipid changes
#   - Pathway Selection: Direct assessment of lipid pathways
#
# Outputs:
#   For each method: hierarchical selection PDFs and Excel result tables
#
################################################################################


################################################################################
# SETUP
## clear workspace and set directories ----
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
################################################################################
# compute run time
run_time <- system.time({
# Set random seed for reproducibility
set.seed(2789465)

# Load required libraries
library(rms)
library(stringr)
library(GSVA)
library(autosgi)
library(tidyverse)
library(magrittr)
library(sgi)
library(maplet)
library(destiny)
library(diffusionMap)

# Load custom functions
source("0_custom_functions.R")


## Define input/output paths ----

# Output directory with timestamp
outdir <- paste0(Sys.Date(), "_adni-results/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Output file prefix
out_prefix <- "adni-lipidomics"


################################################################################
# DATA LOADING
################################################################################

## Load raw ADNI data files ----

# Baker lipidomics data from LONI (longitudinal measurements)
adni_lipidomics_raw <- read.csv("adni-data/ADMCLIPIDOMICSMEIKLELABLONG_08_13_21_20Jun2024.csv") %>% 
  as.data.frame() %>%  
  mutate(uid = paste(RID, VISCODE2, sep = '_'))

# Longitudinal clinical metadata (all timepoints)
clinical_metadata_long <- read_xlsx("adni-data/2024-09-26-long_metadata.xlsx") %>% 
  as.data.frame()

# Baseline clinical metadata
clinical_metadata_base <- read_xlsx("adni-data/2024-09-26-base_metadata.xlsx", sheet = 1) %>% 
  as.data.frame()


################################################################################
# DATA PREPROCESSING
################################################################################

## Clean and format ADNI data ----

# Run comprehensive data cleaning function
# Creates: lipid_se_baseline, colData_wide, clinical_phenotypes
adni_data_cleaning(
  adni_lipidomics_raw, 
  clinical_metadata_long, 
  clinical_metadata_base, 
  outdir, 
  out_prefix
)

# Load cleaned data objects
load(paste0(outdir, out_prefix, "-autosgi-data.rds"))

# Note: The loaded objects are:
#   - se_bl: SummarizedExperiment with baseline lipid data
#   - adni_clinical: Formatted phenotype data for all subjects
#   -  colData_wide:
lipid_se_baseline <- se_bl; clinical_phenotypes <- adni_clinical


################################################################################
# CASE STUDY 1: ssGSEA-BASED AUTOSGI
################################################################################
#
# Approach: Uses single-sample gene set enrichment analysis (ssGSEA) to 
#           compute lipid class-level enrichment scores, then clusters these
#           enrichment patterns to identify groups of lipid classes that 
#           show coordinated changes.
#
################################################################################

## Compute ssGSEA enrichment scores ----

# select lipid classes with more than one lipid annotated to it
lipid_class_groups <- lipid_se_baseline %>% rowData() %>% data.frame() %>% select(class) %>% table() %>% data.frame() %>% 
  filter(Freq>1) %>% pull(class) 

# Calculate pathway enrichment for each lipid class
lipid_ssgsea_scores <- pathway_enrichment(
  dt = (lipid_se_baseline %>% assay() %>% data.frame()), 
  annotations = lipid_class_groups
) %>% 
  as.data.frame() %>% 
  t()

# Assign lipid class names to columns
colnames(lipid_ssgsea_scores) <- lipid_class_groups

# Save enrichment scores for future use
save(lipid_ssgsea_scores, 
     file = paste0(outdir, out_prefix, "-autosgi-ssgsea-data.rds"))


## Run AutoSGI on ssGSEA scores ----

# Ensure rownames match between data and phenotypes
rownames(lipid_ssgsea_scores) <- rownames(clinical_phenotypes)

# Initialize AutoSGI parameters for ssGSEA approach
ssgsea_params <- sgi_params_init(
  dataset = (lipid_ssgsea_scores %>% as.matrix() %>% scale() %>% as.data.frame()), 
  clins = clinical_phenotypes, 
  minsize = (nrow(lipid_ssgsea_scores) / 20), 
  user_defined_tests = c(ordinal = ordinal_test)
)

# Execute hierarchical selection on ssGSEA enrichment scores
ssgsea_results <- hierarchical_selection(
  rule = rule_init(), 
  ssgsea_params, 
  cluster_min = 2, 
  plot = TRUE, 
  supp_plot = TRUE, 
  summary_plot = TRUE, 
  correction_opt = "simes", 
  output_names = list(
    sgi_plots = paste0(outdir, out_prefix, "-ssgsea-hierarchical-selection-results.pdf"), 
    summary = paste0(outdir, out_prefix, "-ssgsea-hierarchical-selection-summary.pdf"), 
    cluster_results = paste0(outdir, out_prefix, "-ssgsea-hierarchical-selection-labels.xlsx"), 
    cluster_indices = paste0(outdir, out_prefix, "-ssgsea-hierarchical-selection-cluster-indices.xlsx"),
    sgi_as_results = paste0(outdir, out_prefix, "-ssgsea-hierarchical-selection-association-results.xlsx")
  )
)
})
print(run_time)
# user  system elapsed 
# 407.755  14.391 401.354 