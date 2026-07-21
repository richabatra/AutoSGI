#Read and preprocess both datasets:

library(rms)
library(stringr)
library(readr)
#get rosmap data
rosmap_data = read_excel("tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx", sheet = "assay")
#get clinical data
rosmap_clins = read_excel("tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx", sheet = "colData")
#get pathways
rosmap_paths = read_excel("tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx", sheet = "rowData")
#transpose rosmap data, correct columns
rosmap_data <- t(rosmap_data)
#use first row (originally header row of assay sheet) as new column names after transpose
colnames(rosmap_data) <- rosmap_data[1, ]
#drop that first row now that it's been used as the header
rosmap_data <- rosmap_data[-1, ]
#define the clinical variables of interest to keep from colData
clin_vars = c("age_death", "msex", "educ", "pmi", "bmi", "anye4", "cogdx", "braaksc", "ceradsc", "sqrt_amyloid", "tangles", "gpath", "niareagansc", "diagnosis", "cogn_global", "cogng_random_slope")
#subset clinical data to only the selected variables
rosmap_clins <- rosmap_clins[, which(colnames(rosmap_clins) %in% clin_vars)]
#coerce to plain data frame
rosmap_clins <- as.data.frame(rosmap_clins)
#convert all assay values to numeric (they were character after transpose)
rosmap_data <- data.frame(apply(rosmap_data, 2, function(x) as.numeric(as.character(x))))
#pull out super pathway annotations for each metabolite
rosmap_super_pathway = rosmap_paths[, 4:4]$SUPER_PATHWAY 
#pull out sub pathway annotations for each metabolite
rosmap_sub_pathway = rosmap_paths[, 5:5]$SUB_PATHWAY 
#align row names of assay data with clinical data row names (subject IDs)
rownames(rosmap_data) <- rownames(rosmap_clins)
#set classes
#APOE e4 carrier status as factor
rosmap_clins$anye4 = as.factor(rosmap_clins$anye4)
#sex as factor
rosmap_clins$msex = as.factor(rosmap_clins$msex)
#NIA-Reagan score as factor
rosmap_clins$niareagansc = as.factor(rosmap_clins$niareagansc)
#diagnosis as factor
rosmap_clins$diagnosis = as.factor(rosmap_clins$diagnosis)
#reverse-code CERAD score so higher values mean more severe pathology (consistent direction with other measures)
rosmap_clins$ceradsc = 5 - rosmap_clins$ceradsc
#mark cognitive diagnosis as ordinal for modeling purposes
class(rosmap_clins$cogdx) = "ordinal"
#mark Braak staging as ordinal
class(rosmap_clins$braaksc) = "ordinal"
#mark (reverse-coded) CERAD score as ordinal
class(rosmap_clins$ceradsc) = "ordinal"
#keep an unmodified copy of the original clinical data for reference
orig_clins = rosmap_clins

#load ADNI lipidomics data
adni_lipids = read.csv("ADMCLIPIDOMICSMEIKLELABLONG_08_13_21_20Jun2024.csv") %>% as.data.frame()

#drop the last column (not needed for analysis)
adni_lipids = adni_lipids[,-ncol(adni_lipids)]
#RIDs are 1057, 1389, 177, 604, 681, 709, 739
#these specific subjects have duplicate m12 visit records that need averaging
multi_rids = c(1057, 1389, 177, 604, 681, 709, 739)
#find rows for these subjects at the m12 visit (the duplicated ones)
midxs = which(adni_lipids$RID %in% multi_rids & adni_lipids$VISCODE2 == "m12")
#remove the duplicate m12 rows, then append one averaged row per subject in their place
adni_lipids = rbind(adni_lipids[-midxs, ], do.call(rbind, lapply(multi_rids, function(rd) {
  #grab all m12 rows for this particular subject
  rows = adni_lipids[which(adni_lipids$RID == rd & adni_lipids$VISCODE2 == "m12"), ] %>% as.data.frame()
  #compute the mean of each lipid measurement column across the duplicate rows
  a = colMeans(rows[,8:ncol(rows)], na.rm = T)
  #build a single collapsed row: keep first row's metadata columns, use averaged lipid values
  dt = data.frame(list(c(as.vector(rows[1:1, 1:7]), colMeans(rows[,8:ncol(rows)], na.rm = T))))
  #restore original column names
  colnames(dt) = colnames(adni_lipids)
  dt
})))
#load long-format ADNI metadata (visit-level clinical outcomes)
colData = read_xlsx("2024-09-26-long_metadata.xlsx") %>% as.data.frame()

#load additional colData (APOE group, education, gender, age) from the q500 metabolomics file
colData2 = read_xlsx("tmp_q500_long.xlsx", sheet = "colData") %>% as.data.frame()

#order tmp_q500 by numbers
#sort colData2 rows by numeric RID so they align properly with colData
o = order(as.numeric(colData2$RID))
colData2 = colData2[o, ]

#we want the patients to have bl outcomes, bl measurement, and be converters
#find subject IDs (RIDs) that have a baseline visit, a non-missing/valid baseline diagnosis (<4), and baseline lipid measurements
patient_rids = Reduce(intersect, list(colData$RID[which(colData$VISCODE == "bl")], colData$RID[which(colData$VISCODE == "bl" & !is.na(colData$DX) & colData$DX < 4)], adni_lipids$RID[which(adni_lipids$VISCODE2 == "bl")]))
#get row indices in colData corresponding to baseline visits for these selected patients
bl_idxs = which(colData$RID %in% patient_rids & colData$VISCODE == "bl")

#subset colData to just these baseline rows to start building the clinical table
adni_clinsl0 = colData[bl_idxs, ]
#add APOE genotype group from colData2
adni_clinsl0$APOEGrp = colData2[bl_idxs, "APOEGrp"]
#add years of education from colData2
adni_clinsl0$PTEDUCAT = colData2[bl_idxs, "PTEDUCAT"]
#add gender from colData2
adni_clinsl0$PTGENDER = colData2[bl_idxs, "PTGENDER"]
#add screening age from colData2
adni_clinsl0$SC_Age = colData2[bl_idxs, "SC_Age"]
#keep only these four demographic/genetic columns for now
adni_clinsl0 = adni_clinsl0[, c("APOEGrp", "PTEDUCAT", "PTGENDER", "SC_Age")]
#helper function to pull all clinical outcome variables for a given visit code (bl, m12, m24) and align them to the baseline subject order
get <- function(id) {
  #initialize an empty (all-NA) placeholder data frame sized to number of patients x 21 outcome columns
  cc = as.data.frame(lapply(seq(21), function(x){rep(NA, times = nrow(adni_clinsl0))}))
  #find which baseline-ordered rows have a record at this visit
  idxs = which(colData$RID[bl_idxs] %in% colData$RID[which(colData$RID %in% patient_rids & colData$VISCODE == id)])
  #extract the actual clinical values for this visit, excluding ID/metadata columns already handled elsewhere
  clin_info = colData[which(colData$RID %in% patient_rids & colData$VISCODE == id), -which(colnames(colData) %in% c("RID", "VISCODE", "BIFAST", "info", "uid", "SC_DXGrp", "APOEGrp", "PTEDUCAT", "PTGENDER", "SC_Age"))]
  #insert these values into the placeholder at the matching subject positions
  cc[idxs, ] = clin_info
  #rename columns to indicate which visit they came from (e.g., "DX_bl", "DX_m12")
  colnames(cc) = sapply(colnames(colData)[c(2, 4:23)], function(x){paste0(x, "_", id)})
  cc
}
#combine baseline demographics with visit-specific outcomes for bl, m12, and m24 visits into one wide table
adni_clinsl0 = cbind(adni_clinsl0, do.call(cbind, lapply(c("bl", "m12", "m24"), get))) %>% as.data.frame()
#convert all columns to numeric (some may have come in as character/factor)
adni_clinsl0 <- data.frame(apply(adni_clinsl0, 2, function(x) as.numeric(as.character(x))))
#mark baseline diagnosis as ordinal
class(adni_clinsl0$DX_bl) = "ordinal"
#mark 12-month diagnosis as ordinal
class(adni_clinsl0$DX_m12) = "ordinal"
#mark 24-month diagnosis as ordinal
class(adni_clinsl0$DX_m24) = "ordinal"
#gender as factor
adni_clinsl0$PTGENDER = as.factor(adni_clinsl0$PTGENDER)
#APOE group as factor
adni_clinsl0$APOEGrp = as.factor(adni_clinsl0$APOEGrp)

#abeta, hip_vol, ptau, tau, adas (sqrt first)
#standardize (z-score) baseline amyloid-beta 42 levels
adni_clinsl0$abeta42_bl <- scale(adni_clinsl0$abeta42_bl) %>% as.numeric()
#standardize 12-month amyloid-beta 42 levels
adni_clinsl0$abeta42_m12 <- scale(adni_clinsl0$abeta42_m12) %>% as.numeric()
#standardize 24-month amyloid-beta 42 levels
adni_clinsl0$abeta42_m24 <- scale(adni_clinsl0$abeta42_m24) %>% as.numeric()
#standardize baseline hippocampal volume
adni_clinsl0$hip_vol_bl <- scale(adni_clinsl0$hip_vol_bl) %>% as.numeric()
#standardize 12-month hippocampal volume
adni_clinsl0$hip_vol_m12 <- scale(adni_clinsl0$hip_vol_m12) %>% as.numeric()
#standardize 24-month hippocampal volume
adni_clinsl0$hip_vol_m24 <- scale(adni_clinsl0$hip_vol_m24) %>% as.numeric()
#standardize baseline phosphorylated tau
adni_clinsl0$ptau_bl <- scale(adni_clinsl0$ptau_bl) %>% as.numeric()
#standardize 12-month phosphorylated tau
adni_clinsl0$ptau_m12 <- scale(adni_clinsl0$ptau_m12) %>% as.numeric()
#standardize 24-month phosphorylated tau
adni_clinsl0$ptau_m24 <- scale(adni_clinsl0$ptau_m24) %>% as.numeric()
#sqrt-transform then standardize baseline ADAS-13 cognitive score (to reduce skew before scaling)
adni_clinsl0$adas13_bl <- scale(adni_clinsl0$adas13_bl %>% sqrt()) %>% as.numeric()
#sqrt-transform then standardize 12-month ADAS-13 score
adni_clinsl0$adas13_m12 <- scale(adni_clinsl0$adas13_m12 %>% sqrt()) %>% as.numeric()
#sqrt-transform then standardize 24-month ADAS-13 score
adni_clinsl0$adas13_m24 <- scale(adni_clinsl0$adas13_m24 %>% sqrt()) %>% as.numeric()

#make sure to first order the lipids by character RID like in clinical measurements (do not have to do this anymore)
#adni_lipids = adni_lipids[order(as.character(adni_lipids$RID)), ]
#subset lipidomics data to only the selected patients at their baseline visit
adni_bl_lipids = adni_lipids[which(adni_lipids$RID %in% patient_rids & adni_lipids$VISCODE2 == "bl"), ]
#drop the first 7 metadata columns, keeping only lipid measurement columns
adni_bl_lipids = adni_bl_lipids[,-seq(1:7)]

#derive lipid class labels by stripping everything after the first "." in each column name
classes = unname(sapply(colnames(adni_bl_lipids), function(c){sub("\\..*", "", c)}))

#function to compute average expression per pathway/class group (simple mean-based aggregation)
pathway_expression <- function(dt, annotations) {
  #get unique group labels
  an = unique(annotations)
  #only keep groups that have more than one member (so averaging makes sense)
  groups = an[sapply(an, function(g){return(length(which(annotations == g)) > 1)})]
  #for each valid group, average across its member rows
  gs = do.call(rbind, lapply(groups, function(g){
    #subset to rows belonging to this group
    subd = dt[which(annotations== g), ]
    #compute column means and transpose into a single row
    t(as.data.frame(colMeans(subd)))
  }))
  gs %>% as.data.frame()
}
#unique lipid class labels
an = unique(classes)
#classes with more than one member (needed for meaningful aggregation)
groups = an[sapply(an, function(g){return(length(which(classes == g)) > 1)})]
#scale lipid data, transpose so lipids are rows, then compute per-class average expression, then transpose back
adni_bl_lipids_avg = as.data.frame(pathway_expression(t(scale(adni_bl_lipids)), classes)) %>% t()
#assign class names as column names
colnames(adni_bl_lipids_avg) = rep(unique(groups), times = 1)

library(GSVA)
#function to compute single-sample GSEA (ssGSEA) enrichment scores per lipid class per sample
pathway_enrichment <- function(dt, annotations) {
  #build gene(lipid) sets: for each class, list the row names (lipid IDs) belonging to it
  gs = lapply(unique(annotations), function(g){return(rownames(dt)[which(annotations== g)])})
  #set up ssGSEA parameters, requiring at least 2 members per set
  param = ssgseaParam(as.matrix(dt), gs, minSize = 2)
  #run ssGSEA
  gsva(param)
}
#transpose so lipids are rows/samples are columns, run ssGSEA enrichment by lipid class, then transpose result back
adni_bl_lipids_ssgsea = as.data.frame(pathway_enrichment(t(adni_bl_lipids), classes)) %>% t()
#unique class labels again
uc = unique(classes)
#restrict to classes with more than one member (matches sets actually used in ssGSEA)
uc = uc[sapply(uc, function(u){length(which(classes == u)) > 1})]
#assign class names as column names of the ssGSEA score matrix
colnames(adni_bl_lipids_ssgsea) = rep(uc, times = 1)

library(ggplot2)
library(cowplot)

#function to generate correlation and hierarchical-clustering-stability plots for a given dataset
make_correlation_plots <- function(
    data,
    dataset_name,
    n_boot = 1000,
    seed = 123
) {
  
  #keep only numeric columns
  data <- data[, vapply(data, is.numeric, logical(1)), drop = FALSE]
  
  #identify columns with nonzero variance (drop constant columns)
  nonzero_variance <- vapply(
    data,
    function(x) sd(x, na.rm = TRUE) > 0,
    logical(1)
  )
  
  #subset to only variable (non-constant) columns
  data <- data[, nonzero_variance, drop = FALSE]
  
  #sanity check: need at least two usable columns to compute correlations
  if (ncol(data) < 2) {
    stop("The dataset must contain at least two numeric, nonconstant columns.")
  }
  
  #z-score standardize all columns
  scaled_data <- scale(data)
  
  #compute pairwise Pearson correlation matrix between all variables
  metabolite_cor_matrix <- cor(
    scaled_data,
    use = "pairwise.complete.obs"
  )
  
  #extract upper triangle (unique pairwise correlations, excluding diagonal)
  metabolite_cors <- metabolite_cor_matrix[upper.tri(metabolite_cor_matrix)]
  #drop any non-finite values (e.g., NA/NaN correlations)
  metabolite_cors <- metabolite_cors[is.finite(metabolite_cors)]
  
  #build a density plot of the pairwise correlation distribution
  correlation_plot <- ggplot(
    data.frame(correlation = metabolite_cors),
    aes(x = correlation)
  ) +
    geom_density(
      fill = "lightblue",
      color = "black",
      alpha = 0.8,
      na.rm = TRUE
    ) +
    labs(
      title = dataset_name,
      x = "Pearson correlation values",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      axis.title.y = element_text(
        angle = 90,
        margin = margin(r = 10)
      ),
      axis.title.x = element_text(
        margin = margin(t = 10)
      )
    )
  
  #fix random seed for reproducible bootstrapping
  set.seed(seed)
  
  #bootstrap procedure: repeatedly resample variables (columns) with replacement and build hierarchical clusterings
  hclusts <- lapply(seq_len(n_boot), function(i) {
    
    #randomly sample column indices with replacement (bootstrap resample of variables)
    sampled_indices <- sample(
      seq_len(ncol(scaled_data)),
      size = ncol(scaled_data),
      replace = TRUE
    )
    
    #build a bootstrap dataset using the resampled columns
    bootstrap_data <- scaled_data[, sampled_indices, drop = FALSE]
    
    #run hierarchical clustering (Ward's method) on this bootstrap sample
    hclust(
      dist(bootstrap_data),
      method = "ward.D2"
    )
  })
  
  #compute cophenetic distances (implied by the dendrogram) for each bootstrap clustering
  cophenetic_distances <- lapply(
    hclusts,
    function(cluster) as.vector(cophenetic(cluster))
  )
  
  #combine all bootstrap cophenetic distance vectors into a matrix (columns = bootstrap replicates)
  cophenetic_matrix <- do.call(cbind, cophenetic_distances)
  
  #compute correlations between cophenetic distance vectors across bootstrap replicates (clustering stability measure)
  hierarchy_cor_matrix <- cor(
    cophenetic_matrix,
    use = "pairwise.complete.obs"
  )
  
  #extract upper triangle of the hierarchy correlation matrix (unique pairwise comparisons)
  hierarchy_cors <- hierarchy_cor_matrix[upper.tri(hierarchy_cor_matrix)]
  #drop non-finite values
  hierarchy_cors <- hierarchy_cors[is.finite(hierarchy_cors)]
  
  #build a density plot of the clustering stability (cophenetic correlation) distribution
  hierarchy_plot <- ggplot(
    data.frame(correlation = hierarchy_cors),
    aes(x = correlation)
  ) +
    geom_density(
      fill = "salmon",
      color = "black",
      alpha = 0.8,
      na.rm = TRUE
    ) +
    labs(
      title = dataset_name,
      x = "Cophenetic distance correlations",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      ),
      axis.title.y = element_text(
        angle = 90,
        margin = margin(r = 10)
      ),
      axis.title.x = element_text(
        margin = margin(t = 10)
      )
    )
  
  #return both plots plus the underlying correlation vectors for further use
  return(
    list(
      correlation_plot = correlation_plot,
      hierarchy_plot = hierarchy_plot,
      metabolite_correlations = metabolite_cors,
      hierarchy_correlations = hierarchy_cors
    )
  )
}
#generate correlation/clustering-stability plots for the ROSMAP metabolomics data
rosmap_plots <- make_correlation_plots(
  data = rosmap_data,
  dataset_name = "ROSMAP",
  n_boot = 1000,
  seed = 123
)

#generate correlation/clustering-stability plots for the ADNI ssGSEA lipid pathway scores
adni_plots <- make_correlation_plots(
  data = as.data.frame(adni_bl_lipids_ssgsea),
  dataset_name = "ADNI",
  n_boot = 1000,
  seed = 123
)
#arrange the four plots (ROSMAP/ADNI correlation and hierarchy plots) into a single 2x2 grid figure
supplementary_figure_4 <- plot_grid(
  rosmap_plots$correlation_plot,
  adni_plots$correlation_plot,
  rosmap_plots$hierarchy_plot,
  adni_plots$hierarchy_plot,
  ncol = 2,
  labels = c("a", "b", "c", "d"),
  align = "hv",
  axis = "tblr"
)

#open a PDF device to save the combined figure
pdf("supplementary_figure_4.pdf",
    width = 12,
    height = 9
)
#draw the figure into the PDF
supplementary_figure_4
#close the PDF device, finalizing the file
dev.off()
