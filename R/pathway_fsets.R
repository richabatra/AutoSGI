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
pathway_fsets <- function(annotations) {
  if (ncol(annotations) != 2) {
    stop("Must have a metabolite colum and an annotation column (1st and 2nd column, respectively)")
  }
  fsets = lapply(unique(annotations[[2]]), function(i){
    if (is.na(i)) {annotations[[1]][which(is.na(annotations[[2]]))]}
    else{annotations[[1]][which(annotations[[2]] == i)]}
  })
  names(fsets) = unique(annotations[[2]])
  fsets
}