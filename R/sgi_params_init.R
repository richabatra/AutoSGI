#' SGI parameter intialization
#'
#' This is a function which initializes the standardized SGI parameters across all AutoSGI runs and checks for data validity
#'
#' @param dataset A dataset with metabolite columns
#' @param clins 
#' @param annotations A dataframe with the first column consisting on the metabolite names and the second column consisting of the corresponding annotations
#' @param method The method to perform reduction, ssgsea or pseudotime
#' @return Always returns an updated dataframe with one column for each annotation
#' @export
#' @examples
#' \dontrun{
#' }
#'
#'
sgi_params_init <- function(dataset, clins, minsize = 18, user_defined_tests = c(), sgi_distance = "euclidean", sgi_linkage = "ward.D2") {
  
  if (nrow(dataset) != nrow(clins)) stop("Number of patient samples in dataset is not equal to number in clins")
  if (!(sgi_distance %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))) stop("Distance must be a valid option from the ones in R dist package")
  if (!(sgi_linkage %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))) stop("Linkage must be a valid option from the ones in R hclust package")
  
  clinical_classes = unique(sapply(clins, class))
  test_classes = c(names(user_defined_tests), "numeric", "factor")
  if (length(setdiff(clinical_classes, test_classes) > 0)) {
    stop("Some clinical variable does not have a supported test (from both user defined tests and default tests)")
  }
  
  if (sum(sapply(dataset, function(x) all(is.na(x) | x == 0))) != 0) {
    warning("A column in dataset is fully empty or constant, this might cause issues with hierarchical clustering")
  }
  if (unique(sapply(dataset, class)) != "numeric") {
    warning("Some columns in dataset are not numeric, this might cause issues with hierarchical clustering")
  }
  if (sum(rownames(dataset) != rownames(clins)) > 0) {
    warning("Rownames of dataset and clinical outcomes not identical, please rename/adjust if possible")
  }
  
  sgi_params <- mget(c("dataset", "clins", "minsize", "user_defined_tests", "sgi_distance", "sgi_linkage"))
  class(sgi_params) <- "sgi_params"
  
  sgi_params
}