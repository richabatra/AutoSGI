set.seed(4956)
# =============================================================================
# Vignette: MOFA+ -> AutoSGI interoperability
#
#   Demonstrates that MOFA2 output plugs into the AutoSGI pipeline with minimal
#   adaptation. MOFA latent factors (and the features behind them) are handed
#   to AutoSGI for subgroup identification and testing, using ROSMAP
#   metabolomics as the worked example.
#
#   MOFA2 is not a competing method here; the point is that its output slots
#   directly into AutoSGI's workflow.
#
# Author: RB
# Date:   2026-07-04
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(MOFA2)
library(tidyverse)
library(magrittr)
library(maplet)
library(autosgi)
library(sgi)
library(ggplot2)

# ---- Custom functions (provides ordinal_test) -------------------------------
source("/case_studies/0_custom_functions.R")

# ---- Input data (ROSMAP, scaled) --------------------------------------------
load('/case_studies/rosmap-results/rosmap-metabo-autosgi-data.rds')

# ---- Output configuration ---------------------------------------------------
# Timestamped output directory and shared file prefix.
outdir <- paste0(Sys.Date(), "_mofa-results/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out_prefix <- "autosgi-mofa"
outfile <- file.path(paste0(outdir, "mofa-model.hdf5"))

# =============================================================================
# 1. Build and train the MOFA model
# =============================================================================

# Create a single-view MOFA object and overwrite its assay with the ROSMAP
# metabolomics matrix.
mobj <- make_example_data(n_views = 1, n_samples = 500, n_features = 667, n_factors = 10)[[1]]
mobj[[1]] <- assay(rosmap_data)
MOFAobject <- create_mofa(mobj)

# Data options.
data_opts <- get_default_data_options(MOFAobject)

# Model options.
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Training options.
train_opts <- get_default_training_options(MOFAobject)

# Prepare and train.
MOFAobject <- prepare_mofa(
  object           = MOFAobject,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)


MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = T)


# =============================================================================
# 2. Downstream analysis and visualisation
# =============================================================================

# Load trained model and attach phenotype metadata.
model <- load_model(outfile)
samples_metadata(model) <- cbind(samples_metadata(model), phenotype_data)

# Variance explained per factor.
plot_variance_explained(model, x = "view", y = "factor")

# Factors 1-3, coloured by cognitive diagnosis.
plot_factor(
  model,
  factors      = c(1, 2, 3),
  color_by     = "cogdx",
  dot_size     = 0.75,      # dot size
  dodge        = TRUE,      # dodge points with different colours
  add_violin   = TRUE,      # add violin plots
  violin_alpha = 0.25       # violin transparency
)


# =============================================================================
# 3. AutoSGI on the MOFA latent factors
# =============================================================================

# Reshape latent factors to a samples x factors matrix.
factors <- get_factors(model, as.data.frame = TRUE)
factors_wide <- factors %>%
  select(sample, factor, value) %>%
  reshape(
    timevar   = "factor",
    idvar     = "sample",
    direction = "wide"
  ) %>%
  data.frame()

# Ensure data and metadata rownames match.
metabolite_matrix <- factors_wide[, -1]
rownames(metabolite_matrix) <- rownames(phenotype_data)

## Hierarchical clustering ----
# Ward.D2 clustering on the latent-space distance matrix.
hierarchical_clusters <- hclust(dist(metabolite_matrix), method = "ward.D2")

## Initialise and run SGI ----
# Minimum cluster size set to 5% of sample size.
sgi_object <- sgi_init(
  hierarchical_clusters,
  minsize            = ceiling(length(hierarchical_clusters$height) / 20),
  outcomes           = phenotype_data,
  user_defined_tests = c(ordinal = ordinal_test)
)

# Execute SGI analysis.
sgi_analysis <- sgi_run(sgi_object)

## SGI visualisation ----
# Tree plot of association results.
sgi_tree_plot <- plot(sgi_analysis, padj_th = 0.05)

# Overview plot with the latent-factor data matrix.
overview_plot <- plot_overview(
  gg_tree = sgi_tree_plot,
  as      = sgi_analysis,
  # outcomes = phenotype_data,  # Commented out per original code
  xdata   = metabolite_matrix
)

# Save plots to PDF.
pdf(
  file   = paste0(outdir, out_prefix, "-mofa-10-latentfactor-sgi.pdf"),
  height = 10,
  width  = 12
)
print(sgi_tree_plot)
print(overview_plot)
dev.off()


# =============================================================================
# 4. Feature sets from the 10 latent factors
# =============================================================================

# Row metadata, keyed by unique feature name.
rd <- rosmap_data %>%
  rowData() %>%
  data.frame() %>%
  mutate(uname = rownames(assay(rosmap_data)))

# For each factor, take the top 10 features by |weight| and label the group.
feat_sets <- list()
for (i in c(1:10)) {
  # Get feature weights, ordered by absolute value.
  W  <- get_weights(model, factors = i, views = 1, as.data.frame = TRUE)
  W1 <- W[order(abs(W$value), decreasing = T), ]
  feat_sets[[i]] <- rd %>%
    filter(uname %in% W1$feature[1:10]) %>%
    mutate(class = paste0("latent_group", i)) %>%
    select(uname, class)
}

# Collapse to one row per feature; retain first factor assignment.
feat_sets <- do.call(rbind, feat_sets)
feat_sets <- feat_sets[-which(duplicated(feat_sets$uname)), ]
feat_sets <- left_join(rd %>% select(uname), feat_sets, by = "uname")


# =============================================================================
# 5. AutoSGI set selection on features behind the latent factors
# =============================================================================

# Raw metabolite matrix, samples x features.
metabolite_matrix <- assay(rosmap_data) %>% t() %>% data.frame()
rownames(metabolite_matrix) <- rownames(phenotype_data)

# Declare ordinal outcomes (natural severity ordering).
class(phenotype_data$cogdx)   <- "ordinal"  # cognitive diagnosis severity
class(phenotype_data$braaksc) <- "ordinal"  # Braak staging (0-VI)
class(phenotype_data$ceradsc) <- "ordinal"  # CERAD score severity

## Hierarchical clustering ----
# Ward.D2 clustering on the latent-space distance matrix.
hierarchical_clusters <- hclust(dist(metabolite_matrix), method = "ward.D2")

## Initialise and run SGI ----
# Minimum cluster size set to 5% of sample size.
sgi_object <- sgi_init(
  hierarchical_clusters,
  minsize            = ceiling(length(hierarchical_clusters$height) / 20),
  outcomes           = phenotype_data,
  user_defined_tests = c(ordinal = ordinal_test)
)

# Execute SGI analysis.
sgi_analysis <- sgi_run(sgi_object)

## SGI visualisation ----
# Tree plot of association results.
sgi_tree_plot <- plot(sgi_analysis, padj_th = 0.05)

# Overview plot with the latent-factor data matrix.
overview_plot <- plot_overview(
  gg_tree = sgi_tree_plot,
  as      = sgi_analysis,
  # outcomes = phenotype_data,  # Commented out per original code
  xdata   = metabolite_matrix
)

# Save plots to PDF.
pdf(
  file   = paste0(outdir, out_prefix, "-mofa-10-latentfactor-sgi.pdf"),
  height = 10,
  width  = 12
)
print(sgi_tree_plot)
print(overview_plot)
dev.off()

## Initialise AutoSGI parameters ----
pathway_params <- sgi_params_init(
  dataset            = metabolite_matrix,
  clins              = phenotype_data,
  minsize            = (nrow(metabolite_matrix) / 20),
  user_defined_tests = c(ordinal = ordinal_test)
)

## Feature-set annotations ----
# Latent-group annotations from the step above.
pathway_annotations <- feat_sets

## Run set selection ----
# Set selection using the latent-factor feature groupings.
pathway_results <- set_selection(
  rule = rule_init(),
  pathway_params,
  pathway_fsets(pathway_annotations),
  plot           = TRUE,
  supp_plot      = TRUE,
  correction_opt = "simes",
  output_names   = list(
    main            = paste0(outdir, out_prefix, "-pathway-selection-results.pdf"),
    cluster_results = paste0(outdir, out_prefix, "-pathway-selection-labels.xlsx"),
    sgi_as_results  = paste0(outdir, out_prefix, "-pathway-selection-association-results.xlsx")
  )
)
