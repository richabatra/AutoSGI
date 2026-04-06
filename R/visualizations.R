#' Pathway outcome plotting
#'
#' @param dataset A dataset with metabolite columns
#' @param clins All clinical outcomes 
#' @param pathways A list of pathway annotations for each metabolite
#' @param pathway The name of the pathway which all the metabolites in the dataset are part of
#' @param by_subgroup Whether to order the plots by subgroup pairs or by phenotypes
#' @param minsize Minimum size for a cluster split to be valid (20 means one twentieth of the entire set)
#' @param sgi_distance Distance metric for representing data before hierarchical clustering
#' @param sgi_linkage Linkage method for hierarchical clustering
#' @param padj_th Subgroup split threshold (and for sgi)
#' @return 1 pdf containing 1. an SGI tree for the pathway and 2. barplots/boxplots for phenotypes
#' @export
pathway_outcomes <- function(dataset, clins, pathways = NULL, pathway = "ALL", by_subgroup = FALSE,
                             minsize = 18, sgi_distance = "euclidean", sgi_linkage = "ward.D2", 
                             padj_th = 0.05) {
  
  #picking the specific pathway from dataset
  if (pathway != "ALL") {
    if (is.null(pathways)) {
      stop("pathways (annotations) cannot be null if pathway is specified")
    }
    pathway_indices = which(pathways == pathway)
    dataset = dataset[, pathway_indices]
  }

  #using set selection to get the objects
  #can just use sgi, but this will run with sgi and we can just extract the plot
  
  #annotations for set selection will just be pathway (essentially all metabolites like sgi)
  cols = length(colnames(dataset))
  annotations = rep(pathway, cols)
  
  #the rule is going to be any hits with specified p value (padj_th)
  rule <- rule_init("ALL", NULL, padj_th)
  #perform set selection
  selection <- set_selection(rule, dataset, clins, annotations, minsize, sgi_distance, sgi_linkage, plot = FALSE, return_plot = TRUE)
  #using set selection so we can get sgi plot
  
  ret_sgi = safe_run_sgi(dataset, sgi_distance, sgi_linkage, minsize, clins)
  
  if (!ret_sgi$invalid & selection$results[1]) { #valid tree
    #run sgi
    sg = ret_sgi$result
    as = sgi_run(sg)
    
    #all significant clusters (ordered by clinical outcome)
    outcomes = colnames(clins)
    sig_clusters = lapply(as$results, extract, as = as, p_value = padj_th)
    
    pdf("phenotype-distributions.pdf")
    par(mar = c(0,0,0,0))
    
    path = "Pathway: "
    path = paste0(path, pathway)
    
    plot.new()
    text(x = 0.5, y = 0.5, path, 
         cex = 0.9, col = "black") 
    
    for (p in selection$plots) {
      plot(p) #plot all sgi objects
      title(main = paste0("Pathway: ", pathway), cex.main = 0.9, line = -1.25)
    }
  
    par(mar = c(0,0,0,0))
    
    summary = "Phenotype Subgroups Distributions\n------------------------------------\n"
    
    max_per_page = 15
    
    if (by_subgroup) {
      
        all_cpairs <- vector()
        for (outcome in outcomes) {
          for (cpair in sig_clusters[outcome]) {
            all_cpairs <- append(all_cpairs, cpair$cluster_pair)
          }
        }
        
        all_cpairs = unique(all_cpairs)
        cpair_results <- vector(length = length(as$cluster_pairs))
        cpair_results <- sapply(rownames(as$cluster_pairs), sig <- function(cpair) {
          return(cpair %in% all_cpairs)
        })
        
        num_subgroups = length(all_cpairs) #number of subgroups having a outcome split
        summary = paste0(summary, num_subgroups, " subgroup pairs out of ", length(rownames(as$cluster_pairs)), " with a significant outcome split\n")
        
        plot(plot_message(summary))
        
        split_data(rownames(as$cluster_pairs), cpair_results, max_per_page, "subgroup-pair")
        
        for (cpair in all_cpairs) {
      
          out = paste0("Subgroup Pair: ", cpair, "\n")
          plot(plot_message(out))
          
          capture.output(plot_outcomes(sg, as, cluster_pairs = c(cpair), padj_th = padj_th))
        }
    }
    else {
      
      sig_cluster <- function(clin, sig_clusters) {
        return(dim(data.frame(sig_clusters[clin]))[1] > 0)
      }
      
      clin_results <- sapply(outcomes, sig_cluster, sig_clusters = sig_clusters)
      
      num_clins = sum(clin_results) #number of clinical outcomes having a subgroup pair
      summary = paste0(summary, num_clins, " clinical outcomes out of ", length(outcomes), " with a significant subgroup pair\n")
      
      plot(plot_message(summary))
      
      #table of phenotypes in subgroup pair?
      split_data(outcomes, clin_results, max_per_page)
      
      for (i in 1:length(clin_results)) {
        if (clin_results[i]) {
          
          out = paste0("Phenotype: ", outcomes[i], "\n")
          plot(plot_message(out))
          
          capture.output(plot_outcomes(sg, as, outcome_names = c(outcomes[i])), padj_th = padj_th) #shouldnt be any nulls
        }
      }
    }
    
    dev.off()
  }
  else {
    pdf("phenotype-distributions.pdf")
    obj = ggplot() + xlim(0, 10) + ylim(0, 10) + theme_void() + 
      annotate("text", x = 5, y = 5, label = "No Results from SGI Analysis", size = 6, colour = "red")
    plot(obj)
    dev.off()
  }
}

#' Correlation outcome plotting
#'
#' @param dataset A dataset with patient rows and clinical outcome columns
#' @return A correlation matrix of the clinical data based on Kendall correlation
#' @export
plot_correlation <- function(clins) {
  
  clin_data <- Filter(function(x) length(unique(x)) > 1, clins)
  
  if (ncol(clin_data) != ncol(clins)) {
    #removing some clinical samples columns -- all same values
    warning("Removing clinical outcomes with patient samples which are all the same")
  }
  if (sum(is.na(clin_data)) != 0) {
    #mean imputation for NAs
    for (i in 1:ncol(clin_data)) {
      clin_data[is.na(clin_data[, i]), i] <- mean(clin_data[, i], na.rm = TRUE)
    }
    warning("Replacing NA values with mean column value")
  }
  
  #kendall correlation and melt
  corr <- cor(clin_data, method = "kendall") %>% 
      as.matrix()
  clin_correlation <- data.frame(Var1 = rownames(corr)[row(corr)], 
                                 Var2=colnames(corr)[col(corr)], 
                                 value=c(corr))
  
  ggplot(data = clin_correlation, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1,1),
                         name="Kendall\nCorrelation\nScale\n") +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
          hjust = 1, size=8), 
          axis.text.y = element_text(size=8)) + 
    xlab("Outcomes") +
    ylab("Outcomes") + 
    ggtitle("Clinical Outcome Correlations") + 
    theme(plot.title = element_text(size = 13, hjust = 0.5, vjust = 1.25)) + 
    coord_fixed() 
  
}