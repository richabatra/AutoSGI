#' Reconstruct an SGI Association Object
#'
#' Reconstructs an SGI association object for a selected cluster from an
#' AutoSGI-generated Excel results file.
#'
#' @param sgi_params The `sgi_params` object used for the original AutoSGI analysis.
#' @param sheet_name Path to the AutoSGI-generated Excel results file.
#' @param cluster_number Number of the cluster to reconstruct.
#' @param cluster_names Character vector of feature names in the cluster.
#' @return An SGI association object containing the adjusted p-values from the results sheet.
#' @export
#' @examples
#' \dontrun{
#' association_object <- get_as(
#'   sgi_params,
#'   "AutoSGI_results.xlsx",
#'   3,
#'   c("feature_1", "feature_2", "feature_3")
#' )
#' }
get_as <- function(sgi_params, sheet_name, cluster_number, cluster_names) {
  as_info = read_xlsx(path = sheet_name, sheet = paste0("SGI ", cluster_number)) %>% as.data.frame()
  rownames(as_info) = as_info[[1]]
  as_info = as_info[,2:ncol(as_info)]
  
  dataset = sgi_params$dataset
  dataset = dataset[,which(colnames(dataset) %in% cluster_names)]
  clins = sgi_params$clins
  minsize = sgi_params$minsize
  user_defined_tests = sgi_params$user_defined_tests
  sgi_distance = sgi_params$sgi_distance 
  sgi_linkage = sgi_params$sgi_linkage
  
  ret_sgi = safe_run_sgi(dataset, minsize, user_defined_tests = user_defined_tests, sgi_distance, sgi_linkage, clins)
  if (class(ret_sgi$result) != "sgi.object") {
    stop("SGI object could not be constructed, invalid cluster")
  }
  
  sg = ret_sgi$result
  as = sgi_run(sg)
  
  lapply(seq(length(as$results)), function(i){
    as$results[[i]][, 3]$padj <<- as_info[[i]]
  })
  
  as
}
