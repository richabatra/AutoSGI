set.seed(7890)
# =============================================================================
# Vignette: ConsensusClusterPlus -> AutoSGI interoperability
#
#   Demonstrates that ConsensusClusterPlus output plugs into the AutoSGI
#   pipeline with minimal adaptation. A consensus clustering tree is coerced
#   into the hclust-shaped object AutoSGI expects, then passed to AutoSGI for
#   subgroup identification and testing. Worked on ROSMAP metabolomics for the
#   full 667-metabolite panel and the 7-metabolite signature, each at k = 2
#   and k = 3.
#
#   ConsensusClusterPlus is not a competing method here; the point is that its
#   output slots directly into AutoSGI's workflow.
#
# Author: RB
# Date:   2026-07-04
# =============================================================================

# ---- One-time install (uncomment if needed) ---------------------------------
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ConsensusClusterPlus")

# ---- Libraries --------------------------------------------------------------
library(ConsensusClusterPlus)
library(tidyverse)
library(magrittr)
library(maplet)
library(autosgi)

# ---- Custom functions (provides ordinal_test) -------------------------------
source("/case_studies/0_custom_functions.R")

# ---- Input data (ROSMAP, scaled) --------------------------------------------
load('/case_studies/rosmap-results/rosmap-metabo-autosgi-data.rds')

## distance method
dist_method <- "pearson" 
# ---- Output paths -----------------------------------------------------------
# Timestamped output directory and shared file prefix.
cc_base <- paste0(Sys.Date(), "_consensus-clustering-results/", dist_method)
dir.create(cc_base, showWarnings = FALSE, recursive = TRUE)
dir_667 <- file.path(cc_base, "667-metabo/")
dir.create(dir_667, showWarnings = FALSE, recursive = TRUE)
dir_7   <- file.path(cc_base, "7-metabo/")
dir.create(dir_7, showWarnings = FALSE, recursive = TRUE)


# ---- Seven-metabolite signature ---------------------------------------------
seven_metabolites <- c(
  "glutamate", "N-acetylglycine", "2-aminoadipate", "guanidinoacetate",
  "glycerophosphoethanolamine", "X - 25020", "glycerophosphorylcholine (GPC)"
)


# =============================================================================
# Helper: extract a consensus tree at a given k and run AutoSGI on it
#   (this is your original per-k block, unchanged, with the varying pieces
#    passed in as arguments)
#
# Note: the tree is always labelled dist.method = "euclidean" / method =
# "ward.D2" to keep the AutoSGI runs comparable to the original case-study
# analyses, regardless of the distance used by ConsensusClusterPlus.
# =============================================================================
run_consensus_sgi <- function(cc_results, k, xdata, minsize_div, outcomes, pdf_file, 
                              dist_method=dist_method) {
  
  # extract hierarchical tree
  hierarchical_clusters             <- cc_results[[k]]$consensusTree
  hierarchical_clusters$labels      <- as.vector(cc_results[[k]]$consensusTree$labels)
  #hierarchical_clusters$dist.method <- as.vector(dist_method)
  hierarchical_clusters$method      <- as.vector("ward.D2")
  
  # Initialize SGI structure
  sgi_object <- sgi_init(
    hierarchical_clusters,
    minsize            = ceiling(length(hierarchical_clusters$height) / minsize_div),
    outcomes           = outcomes,
    user_defined_tests = c(ordinal = ordinal_test)
  )
  
  # Execute SGI analysis
  sgi_analysis <- sgi_run(sgi_object)
  
  ## Generate SGI visualization ----
  sgi_tree_plot <- plot(sgi_analysis, padj_th = 0.2)
  overview_plot <- plot_overview(gg_tree = sgi_tree_plot, as = sgi_analysis, xdata = xdata)
  
  # Save plot to PDF
  pdf(file = pdf_file, height = 8, width = 12)
  print(sgi_tree_plot)
  print(overview_plot)
  dev.off()
  
  invisible(sgi_analysis)
}


# =============================================================================
# 1. Full 667-metabolite panel
# =============================================================================

# feature selection
metabolite_matrix <- assay(rosmap_data) %>% data.frame()
d <- metabolite_matrix %>% as.matrix()

# consensus clustering
results <- ConsensusClusterPlus(
  d, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
  clusterAlg = "hc", distance = dist_method, seed = 1262118388.71279,
  plot = "png", finalLinkage = "ward.D2", title = dir_667
)
save(results, file = file.path(dir_667, "res.rds"))

# minsize divisor per k
#   k=2 -> /20; k=3 -> /50 (the /20 threshold yields no subgroups at k=3,
#   so a smaller minimum size is used to recover results)
minsize_div_667 <- c(`2` = 50, `3` = 50)

# run AutoSGI at k = 2 and k = 3
for (k in 2:3) {
  run_consensus_sgi(
    cc_results  = results,
    k           = k,
    xdata       = t(d),
    minsize_div = minsize_div_667[[as.character(k)]],
    outcomes    = phenotype_data,
    pdf_file    = file.path(dir_667, paste0("cc-667metabo-k", k, "-sgi.pdf"))
  )
}


# =============================================================================
# 2. Seven-metabolite signature
# =============================================================================

# feature selection (subset to the signature)
metabolite_matrix <- metabolite_matrix[which(rowData(rosmap_data)$name %in% seven_metabolites), ]
d <- metabolite_matrix %>% as.matrix()
colnames(d) <- rownames(phenotype_data)

# consensus clustering (Pearson distance; tree still labelled euclidean above
# for comparability with the original analyses)
results <- ConsensusClusterPlus(
  d, maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
  clusterAlg = "hc", distance = dist_method, seed = 1262118388.71279,
  plot = "png", finalLinkage = "ward.D2", title = dir_7
)
save(results, file = file.path(dir_7, "res.rds"))

# minsize divisor per k (both /20)
minsize_div_7 <- c(`2` = 20, `3` = 20)

# run AutoSGI at k = 2 and k = 3
for (k in 2:3) {
  run_consensus_sgi(
    cc_results  = results,
    k           = k,
    xdata       = t(d),
    minsize_div = minsize_div_7[[as.character(k)]],
    outcomes    = phenotype_data,
    pdf_file    = file.path(dir_7, paste0("cc-7metabo-k", k, "-sgi.pdf"))
  )
}
