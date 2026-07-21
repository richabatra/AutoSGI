################################################################################
# ROSMAP METABOLOMICS - AUTOSGI CASE STUDY
################################################################################
#
# Purpose: Analyze ROS/MAP brain metabolomics data using AutoSGI
#          to identify sample clusters associated with clinical phenotypes i.e. subgroups
#
# Study: Religious Orders Study and Memory and Aging Project (ROS/MAP)
# Data: Brain metabolomics corrected for post-mortem interval (PMI)
# Input file is a summarized experiment in an Excel file. This data corresponds to the 
# data from the publication PMID: 35829654. Data is available through sage/synapse. 
# More information on data availablility in the manuscript. 
# The preprocessed data from the above publication was further 
# processed to regress out the effect of post-mortem interval (PMI).
#
# Last modified: Richa Batra
# Date: 2026-07-18
#
# Outputs:
#   - Cleaned data objects (SummarizedExperiment, metabolite matrix, phenotypes)
#   - SGI tree plot showing all associations
#   - AutoSGI hierarchical selection results (PDFs and Excel files)
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
library(autosgi)
library(tidyverse)
library(magrittr)
library(sgi)
library(maplet)

# Load custom functions
source("0_custom_functions.R")


## Define input/output paths ----

# Output directory with timestamp
outdir <- paste0(Sys.Date(), "_rosmap-results/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Output file prefix
out_prefix <- "rosmap-metabo"

# Input data file path
rosmap_data_file <- "rosmap-data/tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx"


################################################################################
# DATA LOADING AND PREPROCESSING
################################################################################

## Load metabolomics data ----

# Load ROS/MAP brain metabolomics data (PMI-corrected)
rosmap_data <- mt_load_se_xls(file = rosmap_data_file)


## Clean and format data ----

# Scale metabolite data for analysis
rosmap_data %<>% mt_pre_trans_scale()

# Extract metabolite matrix (samples × metabolites)
metabolite_matrix <- rosmap_data %>% assay() %>% t() %>% data.frame()

# Extract and format phenotype data
phenotype_data <- rosmap_data %>% colData() %>% data.frame() %>% 
  select(age_death, msex, educ, pmi, bmi, anye4, cogdx, braaksc, ceradsc) %>% 
  
  # Format data types for each variable
  mutate(
    age_death = as.numeric(as.matrix(age_death)),  # Age at death
    msex = as.factor(as.matrix(msex)),             # Sex
    educ = as.numeric(as.matrix(educ)),            # Years of education
    pmi = as.numeric(as.matrix(pmi)),              # Post-mortem interval
    bmi = as.numeric(as.matrix(bmi)),              # Body mass index
    anye4 = as.factor(as.matrix(anye4)),           # APOE4 carrier status
    cogdx = as.factor(as.matrix(cogdx)),           # Cognitive diagnosis
    braaksc = as.factor(as.matrix(braaksc)),       # Braak staging (tau pathology)
    ceradsc = as.factor(as.matrix(ceradsc))        # CERAD score (amyloid pathology)
  )

# Save cleaned data objects for downstream analysis
save(rosmap_data, metabolite_matrix, phenotype_data, 
     file = paste0(outdir, out_prefix, "-autosgi-data.rds"))


################################################################################
# SUBGROUP IDENTIFICATION (SGI) ANALYSIS
################################################################################

## Prepare ordinal variables ----

# Set ordinal class for variables with natural ordering
class(phenotype_data$cogdx) <- "ordinal"    # Cognitive diagnosis severity
class(phenotype_data$braaksc) <- "ordinal"  # Braak staging (0-VI)
class(phenotype_data$ceradsc) <- "ordinal"  # CERAD score severity

# Ensure data and metadata rownames match
rownames(metabolite_matrix) <- rownames(phenotype_data)


## Run hierarchical clustering ----

# Perform hierarchical clustering on metabolite distance matrix
hierarchical_clusters <- hclust(as.dist(metabolite_matrix), method = "ward.D2")


## Initialize and run SGI ----

# Initialize SGI structure
# Minimum cluster size set to 5% of sample size
sgi_object <- sgi_init(
  hierarchical_clusters, 
  minsize = ceiling(length(hierarchical_clusters$height) / 20), 
  outcomes = phenotype_data, 
  user_defined_tests = c(ordinal = ordinal_test)
)

# Execute SGI analysis
sgi_analysis <- sgi_run(sgi_object)


## Generate SGI visualization ----

# Generate tree plot showing association results
sgi_tree_plot <- plot(sgi_analysis, padj_th = 0.5)

# Create overview plot with metabolomics data matrix
overview_plot <- plot_overview(
  gg_tree = sgi_tree_plot, 
  as = sgi_analysis, 
  # outcomes = phenotype_data,  # Commented out per original code
  xdata = metabolite_matrix
)

# Save plot to PDF
pdf(
  file = paste0(outdir, out_prefix, "-all-sgi.pdf"), 
  height = 10, 
  width = 12
)
print(overview_plot)
dev.off()


################################################################################
# AUTOSGI HIERARCHICAL SELECTION ANALYSIS
################################################################################

## Initialize AutoSGI parameters ----

# Set parameters for AutoSGI analysis
autosgi_params <- sgi_params_init(
  dataset = (metabolite_matrix %>% as.data.frame()), 
  clins = phenotype_data, 
  minsize = (nrow(metabolite_matrix) / 20), 
  user_defined_tests = c(ordinal = ordinal_test)
)


## Run hierarchical selection ----

# Execute AutoSGI hierarchical selection procedure
# Generates multiple output files: plots, summaries, and Excel tables
autosgi_results <- hierarchical_selection(
  rule = rule_init(), 
  autosgi_params, 
  cluster_min = 2, 
  plot = TRUE, 
  supp_plot = TRUE, 
  summary_plot = TRUE, 
  correction_opt = "simes", 
  output_names = list(
    sgi_plots = paste0(outdir, out_prefix, "-hierarchical-selection-results.pdf"), 
    summary = paste0(outdir, out_prefix, "-hierarchical-selection-summary.pdf"), 
    cluster_results = paste0(outdir, out_prefix, "-hierarchical-selection-labels.xlsx"), 
    cluster_indices = paste0(outdir, out_prefix, "-hierarchical-selection-cluster-indices.xlsx"),
    sgi_as_results = paste0(outdir, out_prefix, "-hierarchical-selection-association-results.xlsx")
  )
)

})
print(run_time)
# user   system  elapsed 
# 1600.413   48.495 1673.444 
################################################################################
# END OF SCRIPT
################################################################################
