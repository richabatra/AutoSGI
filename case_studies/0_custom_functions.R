################################################################################
# CUSTOM FUNCTIONS FOR AUTOSGI CASE STUDIES
################################################################################
#
# Purpose: This script contains reusable functions for analyzing AutoSGI case studies
#
# Functions included:
#   1. dpt_as_ssgsea()              - Calculate diffusion pseudotime per lipid class
#   2. pathway_enrichment()         - Perform ssGSEA pathway enrichment
#   3. adni_data_cleaning()         - Clean and format ADNI data
#   4. ordinal_test()               - Statistical testing for ordinal outcomes
#   5. get_model_stats()            - Parse ordinal model output
#   6. extract_feature_abundances() - Extract and plot feature abundances
#   7. plot_phenotype_distribution() - Plot phenotype distributions by cluster
#   8. tree_specific_plots()        - Generate tree-specific plots for clusters
#
# Dependencies: DiffusionMap, destiny, ggplot2, dplyr, tidyr, SummarizedExperiment,
#               rms, GSVA, openxlsx, sgi, reshape, reshape2
#
################################################################################


################################################################################
# SECTION 1: DIFFUSION PSEUDOTIME AND PATHWAY ENRICHMENT
################################################################################

#' Calculate Diffusion Pseudotime (DPT) Per Lipid Class
#'
#' This function computes diffusion pseudotime for each lipid class using
#' diffusion maps and diffusion pseudotime analysis.
#'
#' @param eset ExpressionSet object containing lipid data
#' @param grp_list List of lipid classes to analyze
#' @param dist_method Distance metric for diffusion map (default: 'euclidean')
#'
#' @return Matrix of pseudotime values for each lipid class
#'
#' @details Uses a fixed seed (123456) for reproducibility. Applies DiffusionMap
#'          with covertree method for k-nearest neighbors.
#'
dpt_as_ssgsea <- function(eset, grp_list, dist_method = 'euclidean') {
  
  # Loop over each lipid class and calculate pseudotime
  pseudotime_results <- lapply(grp_list, FUN = function(grp) {
    
    # Subset expression set to current lipid class
    current_eset <- eset[which(eset@featureData@data$class %in% grp), ]
    
    # Compute diffusion map with fixed seed for reproducibility
    set.seed(123456)
    diffusion_map <- DiffusionMap(
      current_eset, 
      distance = dist_method, 
      knn_params = list(method = 'covertree')
    )
    
    # Calculate diffusion pseudotime and generate plot
    dpt_object <- diffusion_map %>% DPT() %>% plot()
    
    # Extract pseudotime values from plot object
    pseudotime_values <- dpt_object$data$Colour
    
    return(pseudotime_values)
  }) 
  
  # Combine results across all lipid classes into single matrix
  combined_pseudotime <- do.call(rbind, pseudotime_results)
  
  return(combined_pseudotime)
}


#' Perform ssGSEA Pathway Enrichment Analysis
#'
#' This function performs single-sample Gene Set Enrichment Analysis (ssGSEA)
#' for pathway enrichment using lipid class annotations.
#'
#' @param dt Data matrix with lipids as rows
#' @param annotations Vector of unique lipid class annotations
#'
#' @return ssGSEA enrichment scores
#'
#' @details Creates gene sets from lipid names by extracting prefix before
#'          first period. Minimum gene set size is 2.
#'
pathway_enrichment <- function(dt, annotations) {
  
  # Create annotation dataframe: extract lipid class from rowname prefix
  lipid_annotations <- data.frame(
    name = rownames(dt), 
    annot = sub("\\...*", "", rownames(dt))
  )
  
  # Create gene set list: group lipid names by class
  gene_sets <- lapply(unique(annotations), function(current_class) {
    return(lipid_annotations$name[which(lipid_annotations$annot %in% current_class)])
  })
  names(gene_sets) <- unique(annotations)
  
  # Configure ssGSEA parameters with minimum gene set size of 2
  ssgsea_params <- ssgseaParam(as.matrix(dt), gene_sets, minSize = 2)
  
  # Run ssGSEA enrichment and return scores
  enrichment_scores <- gsva(ssgsea_params)
  
  return(enrichment_scores)
}


################################################################################
# SECTION 2: ADNI DATA CLEANING AND PREPARATION
################################################################################

#' Clean and Format ADNI Data
#'
#' This function performs comprehensive data cleaning for ADNI lipidomics data,
#' including filtering samples, reshaping clinical data, and formatting phenotypes.
#'
#' @param adni_lipids Raw ADNI lipid data
#' @param colData_long Longitudinal clinical metadata (long format)
#' @param colData_base Baseline clinical metadata
#' @param outdir Output directory path
#' @param out_prefix Output file prefix
#'
#' @return Saves cleaned data objects to file: se_bl, colData_wide, adni_clinsl0
#'
#' @details Processes baseline (bl), 12-month (m12), and 24-month (m24) timepoints.
#'          Applies various transformations: scaling for biomarkers, log2 for tau proteins,
#'          sqrt for ADAS13 scores.
#'
adni_data_cleaning <- function(adni_lipids, colData_long, colData_base, outdir, out_prefix) {
  
  ## STEP 1: Filter samples to baseline non-AD convertors ----
  
  # Select samples at baseline with non-missing diagnosis and DX < 4 (non-AD)
  # These are "convertor" subjects who may progress to AD later
  convertors <- colData_long %>% 
    filter(VISCODE %in% "bl") %>% 
    filter(!is.na(DX)) %>%
    mutate(as.numeric(as.matrix(DX))) %>% 
    filter(DX < 4) %>% 
    pull(RID) %>% 
    unique()
  
  
  ## STEP 2: Prepare baseline lipid data ----
  
  # Extract baseline lipid measurements only
  baseline_lipids <- adni_lipids %>% filter(VISCODE2 %in% "bl")
  
  # Set rownames as RID identifiers for easy sample matching
  rownames(baseline_lipids) <- paste0("RID", baseline_lipids$RID)
  
  # Remove non-lipid columns (cols 1-7, 789-790)
  baseline_lipids <- baseline_lipids[, -c(1:7, 789, 790)]
  
  
  ## STEP 3: Create SummarizedExperiment object ----
  
  # Transpose lipids to features (rows) x samples (columns) format
  lipid_matrix <- baseline_lipids %>% t() %>% data.frame()
  
  # Create row metadata: lipid names and classes
  # Class is extracted from lipid name prefix (before first period)
  row_data <- data.frame(
    name = rownames(lipid_matrix), 
    class = sub("\\...*", "", rownames(lipid_matrix))
  )
  
  # Create column metadata: sample IDs joined with clinical data
  col_data <- data.frame(IDs = colnames(lipid_matrix)) %>% 
    mutate(IDs = sub("RID", "", IDs)) %>% 
    mutate(IDs = as.numeric(as.matrix(IDs))) %>%
    left_join(colData_base, by = c("IDs" = "RID")) %>% 
    dplyr::rename(RID = IDs)
  
  # Construct SummarizedExperiment and filter to convertor samples only
  se_bl <- SummarizedExperiment(assay = lipid_matrix, colData = col_data, rowData = row_data) %>%
    mt_modify_filter_samples(RID %in% convertors)
  
  
  ## STEP 4: Identify lipid classes with multiple members ----
  
  # Select lipid classes with more than one lipid (Freq > 1)
  # Single-member classes are excluded from group-level analyses
  lipid_class_groups <- se_bl %>% 
    rowData() %>% 
    data.frame() %>% 
    select(class) %>% 
    table() %>% 
    data.frame() %>% 
    filter(Freq > 1) %>% 
    pull(class)
  
  
  ## STEP 5: Reshape longitudinal clinical data from long to wide format ----
  
  # Reshape clinical data: one row per subject, columns for each timepoint
  # This makes it easier to analyze changes over time
  colData_wide <- reshape(
    colData_long %>% 
      filter(VISCODE %in% c("bl", "m12", "m24")) %>% 
      select(-uid),
    timevar = "VISCODE",
    idvar   = "RID",
    direction = "wide"
  ) %>% 
    data.frame() %>% 
    mutate(RID = as.numeric(as.matrix(RID)))
  
  # Standardize column names (replace "." with "_" for consistency)
  names(colData_wide) <- sub("\\.", "_", names(colData_wide))
  
  
  ## STEP 6: Format phenotype variables ----
  
  # Join clinical data with SummarizedExperiment sample metadata
  colData_wide %<>% 
    right_join((se_bl %>% colData() %>% data.frame() %>% 
                  mutate(RID = as.numeric(as.matrix(RID)))), 
               by = "RID") %>%
    
    ## Format diagnosis and group variables as factors
    mutate(SC_DXGrp_bl = as.factor(as.matrix(SC_DXGrp_bl))) %>%
    mutate(SC_DXGrp_m12 = as.factor(as.matrix(SC_DXGrp_m12))) %>%
    mutate(SC_DXGrp_m24 = as.factor(as.matrix(SC_DXGrp_m24))) %>%
    mutate(DX_bl = as.factor(as.matrix(DX_bl))) %>%
    mutate(DX_m12 = as.factor(as.matrix(DX_m12))) %>%
    mutate(DX_m24 = as.factor(as.matrix(DX_m24))) %>%
    
    ## Format biomarkers: amyloid-beta (scaled)
    mutate(abeta42_bl = as.numeric(scale(as.numeric(as.matrix(abeta42_bl))))) %>%
    mutate(abeta42_m12 = as.numeric(scale(as.numeric(as.matrix(abeta42_m12))))) %>%
    mutate(abeta42_m24 = as.numeric(scale(as.numeric(as.matrix(abeta42_m24))))) %>%
    
    ## Format biomarkers: phospho-tau (log2 transformed and scaled)
    mutate(ptau_bl = as.numeric(scale(log2(as.numeric(as.matrix(ptau_bl)))))) %>%
    mutate(ptau_m12 = as.numeric(scale(log2(as.numeric(as.matrix(ptau_m12)))))) %>%
    mutate(ptau_m24 = as.numeric(scale(log2(as.numeric(as.matrix(ptau_m24)))))) %>%
    
    ## Format biomarkers: total tau (log2 transformed and scaled)
    mutate(tau_bl = as.numeric(scale(log2(as.numeric(as.matrix(tau_bl)))))) %>%
    mutate(tau_m12 = as.numeric(scale(log2(as.numeric(as.matrix(tau_m12)))))) %>%
    mutate(tau_m24 = as.numeric(scale(log2(as.numeric(as.matrix(tau_m24)))))) %>%
    
    ## Format cognitive scores: ADAS13 (sqrt transformed and scaled)
    mutate(adas13_bl = as.numeric(scale(sqrt(as.numeric(as.matrix(adas13_bl)))))) %>%
    mutate(adas13_m12 = as.numeric(scale(sqrt(as.numeric(as.matrix(adas13_m12)))))) %>%
    mutate(adas13_m24 = as.numeric(scale(sqrt(as.numeric(as.matrix(adas13_m24)))))) %>%
    
    ## Format neuroimaging: entorhinal cortex thickness (scaled)
    mutate(ent_th_bl = as.numeric(scale(as.numeric(as.matrix(ent_th_bl))))) %>%
    mutate(ent_th_m12 = as.numeric(scale(as.numeric(as.matrix(ent_th_m12))))) %>%
    mutate(ent_th_m24 = as.numeric(scale(as.numeric(as.matrix(ent_th_m24))))) %>%
    
    ## Format neuroimaging: entorhinal cortex volume (scaled)
    mutate(ent_v_bl = as.numeric(scale(as.numeric(as.matrix(ent_v_bl))))) %>%
    mutate(ent_v_m12 = as.numeric(scale(as.numeric(as.matrix(ent_v_m12))))) %>%
    mutate(ent_v_m24 = as.numeric(scale(as.numeric(as.matrix(ent_v_m24))))) %>%
    
    ## Format neuroimaging: hippocampal volume (scaled)
    mutate(hip_vol_bl = as.numeric(scale(as.numeric(as.matrix(hip_vol_bl))))) %>%
    mutate(hip_vol_m12 = as.numeric(scale(as.numeric(as.matrix(hip_vol_m12))))) %>%
    mutate(hip_vol_m24 = as.numeric(scale(as.numeric(as.matrix(hip_vol_m24))))) %>%
    
    ## Format cognitive composites: ADNI memory score (no scaling)
    mutate(adni_mem_bl = as.numeric(as.matrix(adni_mem_bl))) %>%
    mutate(adni_mem_m12 = as.numeric(as.matrix(adni_mem_m12))) %>%
    mutate(adni_mem_m24 = as.numeric(as.matrix(adni_mem_m24))) %>%
    
    ## Format cognitive composites: ADNI language score (no scaling)
    mutate(adni_lan_bl = as.numeric(as.matrix(adni_lan_bl))) %>%
    mutate(adni_lan_m12 = as.numeric(as.matrix(adni_lan_m12))) %>%
    mutate(adni_lan_m24 = as.numeric(as.matrix(adni_lan_m24))) %>%
    
    ## Format cognitive composites: ADNI executive function score (no scaling)
    mutate(adni_ef_bl = as.numeric(as.matrix(adni_ef_bl))) %>%
    mutate(adni_ef_m12 = as.numeric(as.matrix(adni_ef_m12))) %>%
    mutate(adni_ef_m24 = as.numeric(as.matrix(adni_ef_m24))) %>%
    
    ## Format cognitive composites: ADNI visuospatial score (no scaling)
    mutate(adni_vs_bl = as.numeric(as.matrix(adni_vs_bl))) %>%
    mutate(adni_vs_m12 = as.numeric(as.matrix(adni_vs_m12))) %>%
    mutate(adni_vs_m24 = as.numeric(as.matrix(adni_vs_m24))) %>%
    
    ## Format neuroimaging: global cortical thickness (no scaling)
    mutate(gl_th_bl = as.numeric(as.matrix(gl_th_bl))) %>%
    mutate(gl_th_m12 = as.numeric(as.matrix(gl_th_m12))) %>%
    mutate(gl_th_m24 = as.numeric(as.matrix(gl_th_m24))) %>%
    
    ## Format neuroimaging: global cortical volume (no scaling)
    mutate(gl_v_bl = as.numeric(as.matrix(gl_v_bl))) %>%
    mutate(gl_v_m12 = as.numeric(as.matrix(gl_v_m12))) %>%
    mutate(gl_v_m24 = as.numeric(as.matrix(gl_v_m24))) %>%
    
    ## Format neuroimaging: white matter hyperintensity (no scaling)
    mutate(wmhi_bl = as.numeric(as.matrix(wmhi_bl))) %>%
    mutate(wmhi_m12 = as.numeric(as.matrix(wmhi_m12))) %>%
    mutate(wmhi_m24 = as.numeric(as.matrix(wmhi_m24))) %>%
    
    ## Format neuroimaging: FDG-PET metabolism (no scaling)
    mutate(fdgpet_bl = as.numeric(as.matrix(fdgpet_bl))) %>%
    mutate(fdgpet_m12 = as.numeric(as.matrix(fdgpet_m12))) %>%
    mutate(fdgpet_m24 = as.numeric(as.matrix(fdgpet_m24))) %>%
    
    ## Format clinical measures: BMI (no scaling)
    mutate(SC_BMI_bl = as.numeric(as.matrix(SC_BMI_bl))) %>%
    mutate(SC_BMI_m12 = as.numeric(as.matrix(SC_BMI_m12))) %>%
    mutate(SC_BMI_m24 = as.numeric(as.matrix(SC_BMI_m24))) %>%
    mutate(bmi_bl = as.numeric(as.matrix(bmi_bl))) %>%
    mutate(bmi_m12 = as.numeric(as.matrix(bmi_m12))) %>%
    mutate(bmi_m24 = as.numeric(as.matrix(bmi_m24))) %>%
    
    ## Format neuroimaging: intracranial volume (no scaling)
    mutate(icv_bl = as.numeric(as.matrix(icv_bl))) %>%
    mutate(icv_m12 = as.numeric(as.matrix(icv_m12))) %>%
    mutate(icv_m24 = as.numeric(as.matrix(icv_m24))) %>%
    
    ## Format MRI variables (no scaling)
    mutate(mri_bl = as.numeric(as.matrix(mri_bl))) %>%
    mutate(mri_m12 = as.numeric(as.matrix(mri_m12))) %>%
    mutate(mri_m24 = as.numeric(as.matrix(mri_m24))) %>%
    
    ## Format MRI magnet strength as factor
    mutate(mag_bl = as.factor(as.matrix(mag_bl))) %>%
    mutate(mag_m12 = as.factor(as.matrix(mag_m12))) %>%
    mutate(mag_m24 = as.factor(as.matrix(mag_m24)))
  
  
  ## STEP 7: Format demographic variables ----
  
  colData_wide %<>% 
    mutate(PTGENDER = as.factor(as.matrix(PTGENDER))) %>%
    mutate(APOEGrp = as.factor(as.matrix(APOEGrp))) %>%
    mutate(SC_Age = as.numeric(as.matrix(SC_Age))) %>%
    mutate(PTEDUCAT = as.numeric(as.matrix(PTEDUCAT)))
  
  
  ## STEP 8: Add BMI to SummarizedExperiment and scale lipids ----
  
  # Add baseline BMI as a column in the SummarizedExperiment
  se_bl$bmi_bl <- colData_wide$bmi_bl
  
  # Scale lipid measurements for downstream analysis
  se_bl %<>% mt_pre_trans_scale()
  
  
  ## STEP 9: Create final clinical data object ----
  
  # Select relevant columns for clinical analysis
  # Exclude RID (will be rownames) and redundant diagnostic/MRI columns
  adni_clinical <- colData_wide %>% 
    select(-RID, -SC_DXGrp_bl, -SC_DXGrp_m12, -SC_DXGrp_m24, 
           -mri_bl, -mri_m12, -mri_m24)
  
  # Set rownames as RID for easy sample matching
  rownames(adni_clinical) <- colData_wide$RID
  
  # Set diagnosis variables as ordinal class for proper statistical testing
  class(adni_clinical$DX_bl) <- "ordinal"
  class(adni_clinical$DX_m12) <- "ordinal"
  class(adni_clinical$DX_m24) <- "ordinal"
  
  # Ensure factor variables are properly formatted
  adni_clinical$mag_bl <- as.factor(adni_clinical$mag_bl)
  adni_clinical$mag_m12 <- as.factor(adni_clinical$mag_m12)
  adni_clinical$mag_m24 <- as.factor(adni_clinical$mag_m24)
  adni_clinical$PTGENDER <- as.factor(adni_clinical$PTGENDER)
  adni_clinical$APOEGrp <- as.factor(adni_clinical$APOEGrp)
  
  
  ## STEP 10: Save cleaned data objects ----
  # Save all three cleaned objects to RDS file for downstream analysis
  save(se_bl, colData_wide, adni_clinical, 
       file = paste0(outdir, out_prefix, "-autosgi-data.rds"))
  
  return(invisible(NULL))
}


################################################################################
# SECTION 3: STATISTICAL TESTING FUNCTIONS
################################################################################

#' Ordinal Statistical Test
#'
#' Performs ordinal regression using proportional odds model for ordinal outcomes.
#'
#' @param outcome_var Response variable (ordinal outcome)
#' @param predictor_var Predictor variable
#'
#' @return List containing p-value and test statistic
#'
#' @details Uses rms::orm() with probit family. Returns NA if model fails to converge.
#'          Error handling included via tryCatch.
#'
ordinal_test <- function(outcome_var, predictor_var) {
  
  # Attempt to fit ordinal regression model with error handling
  test_results <- tryCatch({
    
    # Fit proportional odds ordinal regression model
    fitted_model <- rms::orm(outcome_var ~ predictor_var, family = "probit")
    
    # Extract and parse model statistics
    model_stats <- get_model_stats(fitted_model) %>% 
      .$coefs %>% 
      data.frame()
    
    # Standardize column names for clarity
    names(model_stats) <- c('analyte', 'std_error', 'estimate', 
                           'statistic', 'p_value', 'other')
    
    # Return p-value and test statistic for the predictor (last row in output)
    list(
      pval = model_stats$p_value[nrow(model_stats)], 
      stat = model_stats$statistic[nrow(model_stats)]
    )
    
  }, error = function(e) {
    
    # Print error message for debugging purposes
    print(e$message)
    
    # Return NA if test fails (e.g., insufficient samples, convergence issues)
    list(pval = NA, stat = NA)
  })
  
  return(test_results)
}


#' Parse Ordinal Model Output
#'
#' Helper function to extract model statistics from rms ordinal regression output.
#' Code adapted from CoderGuy123 on StackOverflow.
#'
#' @param model_object Fitted ordinal model object from rms::orm()
#' @param precision Number of decimal places to retain (default: 60)
#'
#' @return List containing model statistics and coefficient table
#'
#' @details Temporarily modifies rms::formatNP to preserve full precision,
#'          captures print output, and parses it into structured format.
#'
get_model_stats <- function(model_object, precision = 60) {
  
  ## Preserve full numeric precision in output ----
  
  # Store original formatting function
  old_format_np <- rms::formatNP
  
  # Replace with high-precision formatter temporarily
  assignInNamespace(
    "formatNP", 
    function(x, ...) formatC(x, format = "f", digits = precision), 
    "rms"
  )
  
  # Store original console width setting
  old_width <- options('width')$width
  
  # Increase width to prevent table wrapping during capture
  options(width = old_width + 4 * precision)
  
  
  ## Capture model output as text ----
  
  captured_output <- capture.output(print(model_object))
  
  
  ## Restore original settings ----
  
  options(width = old_width)
  assignInNamespace("formatNP", old_format_np, "rms")
  
  
  ## Extract model statistics ----
  
  model_stats <- c()
  
  # Extract adjusted R-squared if present
  model_stats$R2.adj <- str_match(captured_output, "R2 adj\\s+ (\\d\\.\\d+)") %>% 
    na.omit() %>% 
    .[, 2] %>% 
    as.numeric()
  
  
  ## Extract coefficient table ----
  
  # Find lines containing coefficient information
  coef_start_line <- which(str_detect(captured_output, "Coef\\s+S\\.E\\."))
  coef_end_line <- length(captured_output) - 1
  coef_lines <- captured_output[coef_start_line:coef_end_line]
  
  # Parse coefficient table from captured text
  coef_table <- suppressWarnings(
    readr::read_table(coef_lines %>% stringr::str_c(collapse = "\n"))
  )
  
  # Rename first column for clarity
  colnames(coef_table)[1] <- "Predictor"
  
  # Return structured output
  return(list(stats = model_stats, coefs = coef_table))
}


################################################################################
# SECTION 4: VISUALIZATION FUNCTIONS
################################################################################

#' Extract and Plot Feature Abundances by Cluster
#'
#' Generates boxplots showing metabolite/pathway abundances across clusters.
#'
#' @param method Analysis method: "transform" for transformed data, "pathway" for pathway-level,
#'               or "hc_sel" for hierarchical clustering selection
#' @param cluster_index Cluster index identifier
#' @param main_cl_idx Main cluster index list
#' @param main_df Main data frame
#' @param summarized_exp SummarizedExperiment object with feature data
#' @param path_col Column name for pathway/class annotations (default: "class")
#' @param sample_composition Sample composition across cluster levels
#' @param output_pdf_name Output PDF file path prefix
#'
#' @return Generates PDF with boxplots (no return value)
#'
#' @details Creates separate plots for each cluster pair and pathway/metabolite.
#'          Data is scaled prior to visualization. Uses custom color scheme.
#'
extract_feature_abundances <- function(method = "transform", 
                                      cluster_index, 
                                      main_cl_idx, 
                                      main_df, 
                                      summarized_exp, 
                                      path_col = "class", 
                                      sample_composition,
                                      output_pdf_name) {
  
  ## Define color scheme for boxplots ----
  boxplot_colors <- c("#a8d5ba", "#c3b1e1")
  
  
  ## Scale input data ----
  summarized_exp %<>% mt_pre_trans_scale()
  
  
  ## Prepare data based on analysis method ----
  
  if (method == "transform") {
    
    # Select pathways/features in current cluster
    column_indices <- c(main_cl_idx[[paste0("CLST", cluster_index)]])
    column_indices <- column_indices[which(column_indices > 0)]
    
    # Extract and scale pathway data
    pathway_data <- main_df[, column_indices] %>% 
      as.matrix() %>% 
      scale() %>%
      as.data.frame()
    
    # Store pathway names for plotting
    pathway_names <- names(pathway_data)
    
  } else {
    # For pathway method, use cluster index as pathway name
    pathway_names <- cluster_index
  }
  
  
  ## Initialize PDF output ----
  pdf(
    paste0(output_pdf_name, cluster_index, "_feat_plots.pdf"), 
    width = 3, 
    height = 3
  )
  
  
  ## Loop through each cluster pair ----
  for (cluster_pair in c(setdiff(names(sample_composition), 'sampleids'))) {
    
    # Add title page for current cluster pair
    plot.new()
    text(x = 0.5, y = 0.5, cluster_pair)
    
    # Prepare data with cluster assignments for transformed method
    if (method == "transform") {
      
      cluster_data <- pathway_data %>% 
        data.frame() %>% 
        mutate(clusterids = as.factor(sample_composition[[cluster_pair]]))
      
      # Remove samples without cluster assignment
      cluster_data <- cluster_data %>% 
        .[which(is.na(sample_composition[[cluster_pair]]) == FALSE), ] %>% 
        droplevels()
    }
    
    
    ## Loop through each pathway and create plots ----
    for (pathway in pathway_names) {
      
      # Add title page for current pathway
      plot.new()
      text(x = 0.5, y = 0.5, pathway)
      
      # Plot transformed pathway-level data
      if (method == "transform") {
        
        # Split data by cluster ID for plotting
        pathway_by_cluster <- split(
          cluster_data %>% data.frame() %>% pull(pathway) %>% as.numeric(), 
          cluster_data$clusterids
        )
        
        # Reshape for ggplot2
        plot_data <- reshape::melt(pathway_by_cluster)
        plot_data$value <- as.numeric(plot_data$value)
        plot_data$L1 <- as.factor(plot_data$L1)
        
        # Create boxplot with jittered points overlaid
        pathway_plot <- ggplot(plot_data, aes(y = value, x = L1, fill = L1)) + 
          geom_boxplot() + 
          geom_jitter(alpha = 0.2, size = 0.5) + 
          ylab(pathway) + 
          xlab("cluster ids") + 
          ggtitle(pathway) +
          scale_fill_manual(values = c(boxplot_colors)) +
          theme(legend.title = 'none') +
          theme_bw() + 
          theme(text = element_text(size = 10))
        
        print(pathway_plot)
      }
      
      
      ## Extract metabolites within current pathway ----
      
      if (method == "hc_sel") {
      
        # For hierarchical clustering: select features in this cluster
        column_indices <- c(main_cl_idx[[paste0("CLST", cluster_index)]])
        column_indices <- column_indices[which(column_indices > 0)]
        
        metabolite_data <- main_df[, column_indices] %>% 
          as.matrix() %>% 
          scale() %>%
          as.data.frame()
        
      } else {
        
        # For pathway method: select metabolites annotated to this pathway
        metabolite_indices <- which(
          (summarized_exp %>% rowData() %>% data.frame() %>% pull(path_col)) %in% pathway
        )
        
        metabolite_data <- (summarized_exp %>% assay() %>% t() %>% data.frame()) %>% 
          .[, metabolite_indices] %>%
          as.data.frame()
      }
      
      # Add cluster IDs to metabolite data
      metabolite_cluster_data <- metabolite_data %>% 
        data.frame() %>% 
        mutate(clusterids = as.factor(sample_composition[[cluster_pair]]))
      
      # Remove samples without cluster assignment
      metabolite_cluster_data <- metabolite_cluster_data %>% 
        .[which(is.na(sample_composition[[cluster_pair]]) == FALSE), ] %>% 
        droplevels()
      
      
      ## Plot each individual metabolite ----
      for (metabolite in names(metabolite_data)) {
        
        # Split metabolite values by cluster ID
        metabolite_by_cluster <- split(
          metabolite_cluster_data %>% pull(metabolite) %>% as.numeric(), 
          metabolite_cluster_data$clusterids
        )
        
        # Reshape for ggplot2
        plot_data <- reshape::melt(metabolite_by_cluster)
        plot_data$value <- as.numeric(plot_data$value)
        plot_data$L1 <- as.factor(plot_data$L1)
        
        # Create boxplot with jittered points overlaid
        metabolite_plot <- ggplot(plot_data, aes(y = value, x = L1, fill = L1)) + 
          geom_boxplot() + 
          geom_jitter(alpha = 0.2, size = 0.5) + 
          ylab(metabolite) + 
          xlab("cluster ids") + 
          ggtitle(metabolite) +
          scale_fill_manual(values = c(boxplot_colors)) +
          theme(legend.title = 'none') +
          theme_bw() + 
          theme(text = element_text(size = 10))
        
        print(metabolite_plot)
        
      } # End loop: each metabolite in pathway
      
    } # End loop: each pathway
    
  } # End loop: each cluster pair
  
  # Close PDF device
  dev.off()
  
  return(invisible(NULL))
}


#' Plot Phenotype Distributions Across Clusters
#'
#' Creates visualizations showing how phenotypes are distributed across
#' hierarchical clustering levels.
#'
#' @param output_pdf_name Output PDF file path prefix
#' @param cluster_index Cluster index identifier
#' @param sample_composition Sample composition across cluster levels
#' @param feature_data Feature data frame
#' @param phenotype_data Phenotype data frame
#' @param sig_cluster_pairs Significant cluster pair comparisons
#'
#' @return Generates PDF with phenotype distribution plots (no return value)
#'
#' @details Creates boxplots for numeric phenotypes and stacked bar charts
#'          for categorical/ordinal phenotypes. Only plots significant phenotypes.
#'
plot_phenotype_distribution <- function(output_pdf_name, 
                                       cluster_index, 
                                       sample_composition, 
                                       feature_data,
                                       phenotype_data, 
                                       sig_cluster_pairs) {
  
  ## Define color scheme for plots ----
  categorical_colors <- c("#f08080", "#add8e6", "#d3d3d3", "#d2b48c", 
                         "#ffa07a", "#CCEBC5", "#FBB4AE")
  
  
  ## Initialize PDF output ----
  pdf(
    paste0(output_pdf_name, cluster_index, "_pheno_plots.pdf"), 
    width = 3, 
    height = 3
  )
  
  
  ## Loop over significant clustering levels ----
  for (cluster_pair in c(setdiff(names(sample_composition), 'sampleids'))) {
    
    # Add title page for current cluster level
    plot.new()
    text(x = 0.5, y = 0.5, cluster_pair)
    
    # Subset feature data for samples in this comparison
    cluster_features <- feature_data %>% 
      .[which(is.na(sample_composition[[cluster_pair]]) == FALSE), ]
    
    # Add cluster IDs to phenotype data for grouping
    cluster_phenotypes <- phenotype_data %>% 
      mutate(clusterids = as.factor(sample_composition[[cluster_pair]]))
    
    # Identify phenotypes that are significant at this clustering level
    significant_phenotypes <- sig_cluster_pairs %>% 
      filter(level %in% sub('l', '', cluster_pair)) %>% 
      pull(pheno)
    
    
    ## Loop over significant phenotypes ----
    for (sig_pheno in significant_phenotypes) {
      
      ## Handle numeric phenotypes with boxplots ----
      if (class(cluster_phenotypes[[sig_pheno]]) == 'numeric') {
        
        # Filter to samples with cluster assignments
        pheno_subset <- cluster_phenotypes %>% 
          .[which(is.na(sample_composition[[cluster_pair]]) == FALSE), ] %>% 
          droplevels()
        
        # Split phenotype values by cluster ID
        pheno_by_cluster <- split(
          pheno_subset %>% pull(sig_pheno) %>% as.numeric(), 
          pheno_subset$clusterids
        )
        
        # Reshape for ggplot2
        plot_data <- reshape::melt(pheno_by_cluster)
        plot_data$value <- as.numeric(plot_data$value)
        plot_data$L1 <- as.factor(plot_data$L1)
        
        # Create boxplot showing phenotype distribution by cluster
        numeric_plot <- ggplot(plot_data, aes(y = value, x = L1, fill = L1)) + 
          geom_boxplot() + 
          geom_jitter(alpha = 0.2, size = 0.5) + 
          ylab(sig_pheno) + 
          xlab("cluster ids") + 
          ggtitle(sig_pheno) +
          scale_fill_manual(values = c('#009688', '#F7DC6F')) +
          theme(legend.title = 'none') +
          theme_bw() + 
          theme(text = element_text(size = 10))
        
        print(numeric_plot)
        
      } 
      
      ## Handle categorical/ordinal phenotypes with stacked bars ----
      else if (class(cluster_phenotypes[[sig_pheno]]) %in% c('factor', 'ordinal')) {
        
        # Filter to samples with cluster assignments
        pheno_subset <- cluster_phenotypes %>% 
          .[which(is.na(sample_composition[[cluster_pair]]) == FALSE), ] %>% 
          droplevels()
        
        # Create contingency table: phenotype levels × cluster IDs
        contingency_table <- table(
          pheno_subset %>% pull(sig_pheno) %>% as.numeric(), 
          pheno_subset$clusterids
        )
        
        # Add phenotype level names as rownames for clarity
        rownames(contingency_table) <- pheno_subset %>% 
          pull(sig_pheno) %>% 
          levels()
        
        # Reshape for ggplot2
        plot_data <- reshape::melt(contingency_table)
        plot_data$Var.1 <- as.factor(plot_data$Var.1)
        
        # Create stacked bar chart showing phenotype composition by cluster
        categorical_plot <- ggplot(plot_data, aes(x = Var.2, y = value, fill = Var.1)) + 
          geom_bar(stat = "identity", position = 'fill') +
          ylab(sig_pheno) + 
          xlab("cluster ids") +
          scale_fill_manual(values = categorical_colors) +
          theme(legend.title = 'none') +
          theme_bw() + 
          theme(text = element_text(size = 10))
        
        print(categorical_plot)
        
      } # End conditional: variable type check
      
    } # End loop: significant phenotypes
    
  } # End loop: significant clustering levels
  
  # Close PDF device
  dev.off()
  
  return(invisible(NULL))
}


################################################################################
# SECTION 5: TREE-SPECIFIC ANALYSIS AND PLOTTING
################################################################################

#' Generate Tree-Specific Plots for Clustering Analysis
#'
#' Comprehensive function to run supervised group identification (SGI) analysis
#' and generate multiple types of visualizations for a specific cluster.
#'
#' @param method Analysis method: "transform" or "pathway"
#' @param cluster_of_interest Cluster identifier
#' @param main_df Main data frame with features
#' @param main_cl_idx Main cluster index list
#' @param phenotype_data Phenotype data frame
#' @param padj_threshold Adjusted p-value threshold for significance
#' @param asgi_file Path to existing ASGI results file
#' @param output_pdf_name Output PDF file path prefix
#' @param pathway_column Column name for pathway annotations
#' @param summarized_exp SummarizedExperiment object
#' @param pathway_mapping Optional pathway mapping file (for pathway method)
#'
#' @return Saves SGI results to RDS file and generates PDF plots (no return value)
#'
#' @details Performs hierarchical clustering, runs SGI with custom ordinal test,
#'          compares with existing ASGI results, and generates tree plots and
#'          overview plots.
#'
tree_specific_plots <- function(method,
                               cluster_of_interest, 
                               main_df, 
                               main_cl_idx, 
                               phenotype_data, 
                               padj_threshold, 
                               asgi_file, 
                               output_pdf_name, 
                               pathway_column, 
                               summarized_exp, 
                               pathway_mapping = NULL) {
  
  ## Define color scheme for plots ----
  categorical_colors <- c("#f08080", "#add8e6", "#d3d3d3", "#d2b48c", 
                         "#ffa07a", "#CCEBC5", "#FBB4AE")
  
  
  ## Prepare data based on analysis method ----
  
  if (method == "pathway") {
    
    # Pathway-based analysis approach
    
    # Identify metabolites belonging to current pathway
    metabolite_indices <- which(
      (summarized_exp %>% rowData() %>% 
         data.frame() %>% 
         pull(class)) %in% cluster_of_interest
    )
    
    # Extract metabolite data for this pathway
    cluster_data <- (summarized_exp %>% assay() %>% t() %>% data.frame())[, metabolite_indices]
    
    # Load pathway mappings and corresponding ASGI results
    pathway_map <- read.xlsx(pathway_mapping, sheet = 1)
    asgi_results <- read.xlsx(
      asgi_file, 
      sheet = paste0("SGI ", which(pathway_map$pathways == cluster_of_interest))
    )
    
  } else {
    
    # Transform-based analysis approach
    
    # Select features belonging to current cluster
    column_indices <- c(main_cl_idx[[paste0("CLST", cluster_of_interest)]])
    column_indices <- column_indices[which(column_indices > 0)]
    cluster_data <- main_df[, column_indices]
    
    # Load ASGI results for this cluster
    asgi_results <- read.xlsx(asgi_file, sheet = paste0("SGI ", cluster_of_interest))
  }
  
  
  ## Scale data for clustering ----
  cluster_data <- cluster_data %>% 
    as.matrix() %>% 
    scale() %>%
    as.data.frame()
  
  
  ## Run Supervised Group Identification (SGI) ----
  
  # Ensure data and metadata rownames match
  row.names(cluster_data) <- rownames(phenotype_data)
  
  # Perform hierarchical clustering using Ward's method
  hierarchical_clusters <- hclust(dist(cluster_data), method = "ward.D2")
  
  # Initialize SGI with custom ordinal test
  # Minimum cluster size: 5% of samples (ceiling of height/20)
  sgi_object <- sgi_init(
    hierarchical_clusters, 
    minsize = ceiling(length(hierarchical_clusters$height) / 20), 
    outcomes = phenotype_data, 
    user_defined_tests = c(ordinal = ordinal_test)
  )
  
  # Execute SGI analysis
  sgi_analysis <- sgi::sgi_run(sgi_object)
  
  
  ## Parse and format SGI results ----
  
  # Extract results for all phenotypes and combine into single dataframe
  sgi_results_long <- lapply(1:length(sgi_analysis$results), FUN = function(i) {
    data.frame(sgi_analysis$results[[i]], pheno = names(sgi_analysis$results)[i])
  }) %>% 
    do.call(rbind, .) %>% 
    mutate(CID = paste(cid1, cid2, sep = "vs"))
  
  # Reshape to wide format for comparison with ASGI results
  sgi_results_wide <- sgi_results_long %>% 
    select(CID, pheno, padj) %>% 
    reshape(
      timevar = "pheno",
      idvar   = "CID",
      direction = "wide"
    ) %>% 
    data.frame()
  
  # Clean column names (remove "padj." prefix)
  names(sgi_results_wide) <- sub("padj.", "", names(sgi_results_wide))
  
  
  ## Integrate ASGI significance thresholds ----
  
  # Flag results that meet ASGI significance threshold
  # This allows comparison between current SGI run and previous ASGI results
  for (col_idx in c(1:ncol(sgi_results_wide))) {
    significant_rows <- which(asgi_results[, col_idx] <= padj_threshold)
    sgi_results_wide[significant_rows, col_idx] <- padj_threshold
  }
  
  # Reshape back to long format for easier manipulation
  sgi_results_long_updated <- reshape2::melt(sgi_results_wide) %>% 
    data.frame() %>% 
    dplyr::rename(padj = value) %>%
    mutate(uid = paste0(CID, variable)) %>% 
    select(-CID, -variable)
  
  # Merge with original SGI results to create combined dataset
  sgi_results_combined <- sgi_results_long %<>% 
    mutate(uid = paste0(CID, pheno)) %>% 
    select(-CID) %>% 
    dplyr::rename(padj_sgi = padj) %>%
    left_join(sgi_results_long_updated, by = "uid") %>% 
    select(-uid)
  
  
  ## Identify significant cluster pairs ----
  significant_pairs <- sgi_results_combined %>% filter(padj <= padj_threshold)
  
  
  ## Extract sample cluster assignments ----
  
  # Get virtual cluster paths (VCPs) - shows cluster membership at each hierarchy level
  vcp_data <- get_vcps(sgi_object) %>% data.frame()
  
  # Add sample IDs and select only levels with significant results
  sample_composition <- vcp_data %>% 
    mutate(sampleids = rownames(vcp_data)) %>% 
    select(unique(paste0('l', significant_pairs$level)), sampleids)
  
  
  ## Save results to file ----
  save(
    sample_composition, vcp_data, significant_pairs, 
    sgi_results_long_updated, sgi_analysis, sgi_object, 
    file = paste0(output_pdf_name, cluster_of_interest, "_sgi_data.rds")
  )
  
  
  ## Generate SGI visualization plots ----
  
  pdf(
    paste0(output_pdf_name, cluster_of_interest, "_sgi_plots.pdf"), 
    width = 8, 
    height = 5
  )
  
  # Create SGI tree plot showing hierarchical structure and significance
  sgi_tree_plot <- plot(sgi_analysis, padj_th = padj_threshold)
  
  # Add cluster identifier page
  plot.new()
  text(x = 0.5, y = 0.5, cluster_of_interest)
  
  # Print tree plot
  print(plot(sgi_tree_plot))
  
  # Print overview plot showing tree with phenotypes and feature heatmap
  print(plot_overview(
    gg_tree  = sgi_tree_plot, 
    as       = sgi_analysis, 
    outcomes = phenotype_data, 
    xdata    = cluster_data
  ))
  
  dev.off()
  
  
  ## Generate feature abundance plots ----
  extract_feature_abundances(
    method, cluster_of_interest, main_cl_idx, main_df, summarized_exp, 
    pathway_column, sample_composition, output_pdf_name
  )
  
  
  ## Generate phenotype distribution plots ----
  plot_phenotype_distribution(
    output_pdf_name, cluster_of_interest, sample_composition, 
    cluster_data, phenotype_data, significant_pairs
  )
  
  return(invisible(NULL))
}


################################################################################
# END OF SCRIPT
################################################################################
