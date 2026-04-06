#' Pathway Feature Set Construction
#'
#' This is a function which constructs feature sets for AutoSGI based on biomolecular pathway annotations
#'
#' @param annotations A dataframe with the first column consisting on the metabolite names and the second column consisting of the corresponding annotations
#' @return Returns a list of lists with a set of metabolite names in each list
#' @export
#' @examples
#' \dontrun{
#' }
#'
#'
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