#' Set selection
#'
#' This is a function which uses predefined features sets of metabolites and finds \code{sgi} trees of interest
#' based on a \code{rule} object
#'
#' @param rule A \code{rule} object containing a target clinical outcome, maximum level, and p value
#' @param sgi_params A \code{sgi_params} object containing information including the data, clinical outcomes, etc
#' @param fsets A list of lists containing the feature sets to run AutoSGI on. The names of the sets will be inferred from the names of the lists. 
#' @param correction_opt The method used for multiple testing corrected, either simes or bonferroni
#' @param plot Whether to display plots of sgi trees that are of interest according to the selection rule
#' @param supp_plot Whether to display all the feature names of the pathway cluster in the plot file
#' @param output_names The output names of the AutoSGI files
#' @return A vector with a boolean value for each corresponding metabolite pathway, indicating whether the sgi tree formed by that subset of metabolites is of interest
#' @export
#' @examples
#' \dontrun{
#' dataset <- autosgi::qmdiab_dataset
#' annotations <- autosgi::qmdiab_annotations
#' #two categorizations of pathways
#' super_pathway <- annotations[, 3:3]$SUPER_PATHWAY #subsetting on metabolite rows
#' sub_pathway <- annotations[, 4:4]$SUB_PATHWAY
#' rule <- rule_init(target_outcome = "BMI", max_level = 3, p_value = 0.05) #random rule, can change based on dataset
#' two types of pathway annotations for this case
#' sup_tests <- set_selection(rule, sgi::qmdiab_plasma, sgi::qmdiab_clin, super_pathway, minsize = 18, plot = FALSE)
#' sub_tests <- set_selection(rule, sgi::qmdiab_plasma, sgi::qmdiab_clin, sub_pathway, minsize = 18, plot = FALSE)
#' print(sup_tests)
#' print(sub_tests)
#' }
#'
set_selection <- function(rule, sgi_params, fsets, correction_opt = "simes", plot = T, supp_plot = F,
                          output_names = list(main = "set-selection-results.pdf",
                                              cluster_results = "set-selection-labels.xlsx",
                                              sgi_as_results = "set-selection-association-results.xlsx")) {
  
  #get information from sgi_params object
  dataset = sgi_params$dataset
  clins = sgi_params$clins
  minsize = sgi_params$minsize
  user_defined_tests = sgi_params$user_defined_tests
  sgi_distance = sgi_params$sgi_distance 
  sgi_linkage = sgi_params$sgi_linkage
  
  #check data
  if (correction_opt != "simes" & correction_opt != "bonferroni" & correction_opt != "extended_simes") stop("correction must be extended simes, simes, or bonferroni.")
  fnames = unique(unlist(fsets))
  if (length(setdiff(fnames, colnames(dataset)) > 0)) {
    stop("Some feature names are not part of dataset column names")
  }
  
  if (is.null(names(fsets))) {
    names(fsets) = paste0("Cluster-", seq(length(fsets)))
  }

  #get information from rule object
  target_outcome = rule$target_outcome
  max_level = rule$max_level
  p_value = rule$p_value
  
  #filter out fsets with only one metabolite
  #fsets = fsets[seq(length(fsets))[lapply(fsets, length) > 1]]
  if (min(sapply(fsets, length)) == 1) {
    print("Removing feature sets with length 1")
  }
  unique_annotations = names(fsets)
  unique_annotations[is.na(unique_annotations)] = "NA"
  
  #go through all annotation categories
  plots = vector(mode = "list")
  
  #errors wont matter until we get to plotting (unless they affect all the sgis)
  #collect sgs first
  sgs <- lapply(fsets, function(fset) {
    idxs = which(colnames(dataset) %in% fset)
    metabolite_cluster <- dataset[, idxs] %>% as.data.frame()
    ret_sgi = safe_run_sgi(metabolite_cluster, minsize, user_defined_tests = user_defined_tests, sgi_distance, sgi_linkage, clins)
    ret_sgi$result
  })
  #run association testing
  asx <- lapply(sgs, function(sg) {
    tryCatch({
      sgi_run(sg)
    },
    error = function(e){
      print(e$message)
    })
  })
  #correct everything
  asx <- correct_as(sgs, asx, clins, correction_opt)
  #now we check validation
  res <- lapply(seq(length(unique_annotations)), function(i) {
    #some var initializing
    inv_clust <- val_tree <- F 
    annotation = unique_annotations[i]
    
    if (tree_validation(asx[[i]], rule)) {
      val_tree <- T
    }
    
    if (plot) {
      
      if (nchar(annotation) > 6) {
        annotation <- paste("path", as.character(i))
      }
      
      #run all methods needed for plotting
      #if theres an error in sgi, cluster is invalid
      #use run_plot (try-catch) for this
      ret_plot = safe_run_plot(asx[[i]], p_value, clins, dataset[, which(colnames(dataset) %in% fsets[[i]])], annotation)
      plots[[length(plots) + 1]] <<- ret_plot$result
      inv_clust <- ret_plot$invalid
    }
    
    if (inv_clust) {
      label <- "Invalid"
      val_tree <- F
    }
    else {
      if (val_tree) {label <- "True"}
      else {label <- "False"}
    }
    return(list(valid_tree = val_tree, label = label))
  })
  res <- as.data.frame(do.call(rbind, res))
  valid_tree = av(res$valid_tree)
  labels = av(res$label)
  
  if (plot) {
    
    pdf(output_names$main, onefile = TRUE)
    par(mar = c(0,0,0,0))
    summary = "Set-Selection Plots\n------------------------------------\n"
    num_interest = sum(valid_tree) #number of trees satisfying rule
    summary = paste0(summary, num_interest, " pathways of interest, out of ", length(unique_annotations), " pathways.\n")
    summary = paste0(summary, "Rule used: ", capture.output(print.rule(rule)), "\n")
    plot.new()
    
    text(x = 0.5, y = 0.5, paste(summary), 
         cex = 0.9, col = "black")
    
    max_per_page = 15
    split_data(unique_annotations, valid_tree, max_per_page)
    
    par(oma=c(0, 0, 0, 0))
    par(mar=c(1.1, 1.1, 1.1, 1.1))
    
    sapply(seq(length(unique_annotations)), function(i){
      plot(plots[[i]])
      title(main = paste0("Pathway: ", unique_annotations[i]), cex.main = 0.9, line = -0.5)
      #supp_plot with tree null 
      if (supp_plot) {
        display_features(fsets[[i]], maximum = 10)
      }
    })
    
    dev.off()
    
    #also return an excel file with as info
    cluster_data <- data.frame(pathways = unique_annotations, labels = labels)
    write_xlsx(cluster_data, output_names$cluster_results)
    
    as_sheet = combined_as_sheet(asx)
    names(as_sheet) = do.call("c", lapply(seq(length(asx)), function(i){paste0("SGI ", i)}))
    write_xlsx(as_sheet, path = output_names$sgi_as_results)
  }
  
  return(list(results = valid_tree))
}
