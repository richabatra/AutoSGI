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
colnames(rosmap_data) <- rosmap_data[1, ]
rosmap_data <- rosmap_data[-1, ]
clin_vars = c("age_death", "msex", "educ", "pmi", "bmi", "anye4", "cogdx", "braaksc", "ceradsc", "sqrt_amyloid", "tangles", "gpath", "niareagansc", "diagnosis", "cogn_global", "cogng_random_slope")
rosmap_clins <- rosmap_clins[, which(colnames(rosmap_clins) %in% clin_vars)]
rosmap_clins <- as.data.frame(rosmap_clins)
rosmap_data <- data.frame(apply(rosmap_data, 2, function(x) as.numeric(as.character(x))))
rosmap_super_pathway = rosmap_paths[, 4:4]$SUPER_PATHWAY 
rosmap_sub_pathway = rosmap_paths[, 5:5]$SUB_PATHWAY 
rownames(rosmap_data) <- rownames(rosmap_clins)
#set classes
rosmap_clins$anye4 = as.factor(rosmap_clins$anye4)
rosmap_clins$msex = as.factor(rosmap_clins$msex)
rosmap_clins$niareagansc = as.factor(rosmap_clins$niareagansc)
rosmap_clins$diagnosis = as.factor(rosmap_clins$diagnosis)
rosmap_clins$ceradsc = 5 - rosmap_clins$ceradsc
class(rosmap_clins$cogdx) = "ordinal"
class(rosmap_clins$braaksc) = "ordinal"
class(rosmap_clins$ceradsc) = "ordinal"
orig_clins = rosmap_clins

adni_lipids = read.csv("ADMCLIPIDOMICSMEIKLELABLONG_08_13_21_20Jun2024.csv") %>% as.data.frame()

adni_lipids = adni_lipids[,-ncol(adni_lipids)]
#RIDs are 1057, 1389, 177, 604, 681, 709, 739
multi_rids = c(1057, 1389, 177, 604, 681, 709, 739)
midxs = which(adni_lipids$RID %in% multi_rids & adni_lipids$VISCODE2 == "m12")
adni_lipids = rbind(adni_lipids[-midxs, ], do.call(rbind, lapply(multi_rids, function(rd) {
  rows = adni_lipids[which(adni_lipids$RID == rd & adni_lipids$VISCODE2 == "m12"), ] %>% as.data.frame()
  a = colMeans(rows[,8:ncol(rows)], na.rm = T)
  dt = data.frame(list(c(as.vector(rows[1:1, 1:7]), colMeans(rows[,8:ncol(rows)], na.rm = T))))
  colnames(dt) = colnames(adni_lipids)
  dt
})))
colData = read_xlsx("2024-09-26-long_metadata.xlsx") %>% as.data.frame()

colData2 = read_xlsx("tmp_q500_long.xlsx", sheet = "colData") %>% as.data.frame()

#order tmp_q500 by numbers
o = order(as.numeric(colData2$RID))
colData2 = colData2[o, ]

#we want the patients to have bl outcomes, bl measurement, and be converters
patient_rids = Reduce(intersect, list(colData$RID[which(colData$VISCODE == "bl")], colData$RID[which(colData$VISCODE == "bl" & !is.na(colData$DX) & colData$DX < 4)], adni_lipids$RID[which(adni_lipids$VISCODE2 == "bl")]))
bl_idxs = which(colData$RID %in% patient_rids & colData$VISCODE == "bl")

adni_clinsl0 = colData[bl_idxs, ]
adni_clinsl0$APOEGrp = colData2[bl_idxs, "APOEGrp"]
adni_clinsl0$PTEDUCAT = colData2[bl_idxs, "PTEDUCAT"]
adni_clinsl0$PTGENDER = colData2[bl_idxs, "PTGENDER"]
adni_clinsl0$SC_Age = colData2[bl_idxs, "SC_Age"]
adni_clinsl0 = adni_clinsl0[, c("APOEGrp", "PTEDUCAT", "PTGENDER", "SC_Age")]
get <- function(id) {
  cc = as.data.frame(lapply(seq(21), function(x){rep(NA, times = nrow(adni_clinsl0))}))
  idxs = which(colData$RID[bl_idxs] %in% colData$RID[which(colData$RID %in% patient_rids & colData$VISCODE == id)])
  clin_info = colData[which(colData$RID %in% patient_rids & colData$VISCODE == id), -which(colnames(colData) %in% c("RID", "VISCODE", "BIFAST", "info", "uid", "SC_DXGrp", "APOEGrp", "PTEDUCAT", "PTGENDER", "SC_Age"))]
  cc[idxs, ] = clin_info
  colnames(cc) = sapply(colnames(colData)[c(2, 4:23)], function(x){paste0(x, "_", id)})
  cc
}
adni_clinsl0 = cbind(adni_clinsl0, do.call(cbind, lapply(c("bl", "m12", "m24"), get))) %>% as.data.frame()
adni_clinsl0 <- data.frame(apply(adni_clinsl0, 2, function(x) as.numeric(as.character(x))))
class(adni_clinsl0$DX_bl) = "ordinal"
class(adni_clinsl0$DX_m12) = "ordinal"
class(adni_clinsl0$DX_m24) = "ordinal"
adni_clinsl0$PTGENDER = as.factor(adni_clinsl0$PTGENDER)
adni_clinsl0$APOEGrp = as.factor(adni_clinsl0$APOEGrp)

#abeta, hip_vol, ptau, tau, adas (sqrt first)
adni_clinsl0$abeta42_bl <- scale(adni_clinsl0$abeta42_bl) %>% as.numeric()
adni_clinsl0$abeta42_m12 <- scale(adni_clinsl0$abeta42_m12) %>% as.numeric()
adni_clinsl0$abeta42_m24 <- scale(adni_clinsl0$abeta42_m24) %>% as.numeric()
adni_clinsl0$hip_vol_bl <- scale(adni_clinsl0$hip_vol_bl) %>% as.numeric()
adni_clinsl0$hip_vol_m12 <- scale(adni_clinsl0$hip_vol_m12) %>% as.numeric()
adni_clinsl0$hip_vol_m24 <- scale(adni_clinsl0$hip_vol_m24) %>% as.numeric()
adni_clinsl0$ptau_bl <- scale(adni_clinsl0$ptau_bl) %>% as.numeric()
adni_clinsl0$ptau_m12 <- scale(adni_clinsl0$ptau_m12) %>% as.numeric()
adni_clinsl0$ptau_m24 <- scale(adni_clinsl0$ptau_m24) %>% as.numeric()
adni_clinsl0$adas13_bl <- scale(adni_clinsl0$adas13_bl %>% sqrt()) %>% as.numeric()
adni_clinsl0$adas13_m12 <- scale(adni_clinsl0$adas13_m12 %>% sqrt()) %>% as.numeric()
adni_clinsl0$adas13_m24 <- scale(adni_clinsl0$adas13_m24 %>% sqrt()) %>% as.numeric()

#make sure to first order the lipids by character RID like in clinical measurements (do not have to do this anymore)
#adni_lipids = adni_lipids[order(as.character(adni_lipids$RID)), ]
adni_bl_lipids = adni_lipids[which(adni_lipids$RID %in% patient_rids & adni_lipids$VISCODE2 == "bl"), ]
adni_bl_lipids = adni_bl_lipids[,-seq(1:7)]

classes = unname(sapply(colnames(adni_bl_lipids), function(c){sub("\\..*", "", c)}))

pathway_expression <- function(dt, annotations) {
  an = unique(annotations)
  groups = an[sapply(an, function(g){return(length(which(annotations == g)) > 1)})]
  gs = do.call(rbind, lapply(groups, function(g){
    subd = dt[which(annotations== g), ]
    t(as.data.frame(colMeans(subd)))
  }))
  gs %>% as.data.frame()
}
an = unique(classes)
groups = an[sapply(an, function(g){return(length(which(classes == g)) > 1)})]
adni_bl_lipids_avg = as.data.frame(pathway_expression(t(scale(adni_bl_lipids)), classes)) %>% t()
colnames(adni_bl_lipids_avg) = rep(unique(groups), times = 1)

library(GSVA)
pathway_enrichment <- function(dt, annotations) {
  gs = lapply(unique(annotations), function(g){return(rownames(dt)[which(annotations== g)])})
  param = ssgseaParam(as.matrix(dt), gs, minSize = 2)
  gsva(param)
}
adni_bl_lipids_ssgsea = as.data.frame(pathway_enrichment(t(adni_bl_lipids), classes)) %>% t()
uc = unique(classes)
uc = uc[sapply(uc, function(u){length(which(classes == u)) > 1})]
colnames(adni_bl_lipids_ssgsea) = rep(uc, times = 1)

library(ggplot2)
library(cowplot)

make_correlation_plots <- function(
    data,
    dataset_name,
    n_boot = 1000,
    seed = 123
) {
  
  data <- data[, vapply(data, is.numeric, logical(1)), drop = FALSE]
  
  nonzero_variance <- vapply(
    data,
    function(x) sd(x, na.rm = TRUE) > 0,
    logical(1)
  )
  
  data <- data[, nonzero_variance, drop = FALSE]
  
  if (ncol(data) < 2) {
    stop("The dataset must contain at least two numeric, nonconstant columns.")
  }
  
  scaled_data <- scale(data)
  
  metabolite_cor_matrix <- cor(
    scaled_data,
    use = "pairwise.complete.obs"
  )
  
  metabolite_cors <- metabolite_cor_matrix[upper.tri(metabolite_cor_matrix)]
  metabolite_cors <- metabolite_cors[is.finite(metabolite_cors)]
  
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
  
  set.seed(seed)
  
  hclusts <- lapply(seq_len(n_boot), function(i) {
    
    sampled_indices <- sample(
      seq_len(ncol(scaled_data)),
      size = ncol(scaled_data),
      replace = TRUE
    )
    
    bootstrap_data <- scaled_data[, sampled_indices, drop = FALSE]
    
    hclust(
      dist(bootstrap_data),
      method = "ward.D2"
    )
  })
  
  cophenetic_distances <- lapply(
    hclusts,
    function(cluster) as.vector(cophenetic(cluster))
  )
  
  cophenetic_matrix <- do.call(cbind, cophenetic_distances)
  
  hierarchy_cor_matrix <- cor(
    cophenetic_matrix,
    use = "pairwise.complete.obs"
  )
  
  hierarchy_cors <- hierarchy_cor_matrix[upper.tri(hierarchy_cor_matrix)]
  hierarchy_cors <- hierarchy_cors[is.finite(hierarchy_cors)]
  
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
  
  return(
    list(
      correlation_plot = correlation_plot,
      hierarchy_plot = hierarchy_plot,
      metabolite_correlations = metabolite_cors,
      hierarchy_correlations = hierarchy_cors
    )
  )
}
rosmap_plots <- make_correlation_plots(
  data = rosmap_data,
  dataset_name = "ROSMAP",
  n_boot = 1000,
  seed = 123
)

adni_plots <- make_correlation_plots(
  data = as.data.frame(adni_bl_lipids_ssgsea),
  dataset_name = "ADNI",
  n_boot = 1000,
  seed = 123
)
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

pdf("supplementary_figure_4.pdf",
    width = 12,
    height = 9
)
supplementary_figure_4
dev.off()