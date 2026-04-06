#' Validate a \code{sgi} tree based on selection rule
#'
#' This is a function which determines if a tree is of interest based on a \code{rule} object
#'
#' @param sgi_tree an \code{as} object generated from metabolite features and clinical outcomes
#' @param target_outcome the clinical outcome that is of interest
#' @param max_level the maximum depth that the first target outcome can have to be of interest
#' @param p_value the p value used for association testing
#' @return A boolean value, with \code{TRUE} denoting a tree of interest and \code{FALSE} otherwise, all based on the selection rule
#' @export
#' @examples
#' \dontrun{
#' hc = hclust(dist(sgi::qmdiab_plasma), method = "ward.D2")
#' sg = sgi_init(hc, minsize = 18, outcomes = sgi::qmdiab_clin)
#' as = sgi_run(sg)
#' rule1 <- rule_init("BMI", 3, 0.05)
#' rule2 <- rule_init("BIM", 3, 0.05)
#' rule3 <- rule_init("Cholesterol", 3, 0.05)
#' rule4 <- rule_init("Cholesterol", 3, 0.01)
#' rule5 <- rule_init("Diabetes", 2, 0.05)
#' tree_validation(as, rule1)
#' tree_validation(as, rule2)
#' tree_validation(as, rule3)
#' tree_validation(as, rule4)
#' tree_validation(as, rule5)
#' }
#'
tree_validation <- function(as, rule) {
  
  #default targeted_outcome is all clinical outcomes
  #default max_level is any tree depth
  #default p_value (padj_th) is 0.05
  
  #get information from rule object
  target_outcome = rule$target_outcome
  max_level = rule$max_level
  p_value = rule$p_value
  
  valid_tree = FALSE #stores output (whether a tree is valid based on the selection rule)
  
  full_clins <- names(as$results)
  clinical_outcomes <- full_clins
  
  if (target_outcome != "ALL") {
    #use all possible clinical outcomes
    clinical_outcomes <- c(target_outcome)
  }
  
  #all outcome clusters
  outcome_clusters = lapply(as$results, extract, as = as, p_value = p_value)
  
  sapply(clinical_outcomes, function(outcome) {
    data <- data.frame(outcome_clusters[outcome])
    #check if a targeted outcome is in the entire clinical outcomes
    if (dim(data)[1] == 0) {
      if (outcome %in% full_clins) {
        return()
      }
      else {
        stop("one of targeted outcomes is not in clinical outcomes")
      }
    }
    clusters <- data[, 1:1]
    levels <- data[, 6:6]
    sapply(seq(length(clusters)), function(index){
      l = levels[index]
      
      #multiple checks for validity
      #l should not be NA
      #l should be at most max_level for a valid tree
      #if max_level is null l just has to be valid
      
      if (l == 'NA') {
        return()
      }
      c = clusters[index]
      
      #within max_level depth
      if (isTRUE(l <= max_level)) {
        valid_tree <<- T
      }
      else if (is.null(max_level)) { #max_level is any depth
        valid_tree <<- T
      }
    })
  })

  valid_tree #return output
  
}
