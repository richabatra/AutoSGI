#'Initialize \code{rule} object
#'
#'This function initializes a \code{rule} object to be used to check if an \code{sgi} tree is of interest
#'
#' @param target_outcome the clinical outcome that is of interest for "splitting" a cluster into two subgroups
#' @param max_level the maximum depth that the first clinical_outcome split can have to be of interest
#' @param p_value p value for association testing
#' @return a \code{rule} object which can be used for tree validation
#' @export
rule_init <- function(target_outcome = "ALL", max_level = NULL, p_value = 0.05) {
  
  #rule.init adds the target outcome, maximum level, and p value to a rule object
  
  if (p_value < 0 | p_value > 1) {
    stop("p value has to be in the range [0, 1]")
  }
  
  if (!is.null(max_level)) {
    if (max_level <= 1) {
      stop("max_level must be at least 2")
    }
  }
  
  rule <- list(target_outcome = target_outcome, max_level = max_level, p_value = p_value) #creating object
  class(rule) <- "rule"
  
  rule
}

#'Print \code{rule} object summary
#'
#'This function prints a summary of all the variables contained within a rule object
#'
#' @param rule an object of the class \code{rule}
#' @return a description of the clinical outcome, max level, and p value contained in the object
#' @export
#' 
#' @examples
#' \dontrun{
#' rule1 <- rule_init("BMI", 3, 0.05)
#' rule2 <- rule_init("BIM", 3, 0.1)
#' rule3 <- rule_init("Cholesterol", 3, -1) #invalid
#' rule4 <- rule_init("Cholesterol", 3, 1.00001) #invalid
#' rule5 <- rule_init("Diabetes", 2, 0.05)
#' print(rule1)
#' print(rule2)
#' print(rule3)
#' print(rule4)
#' print(rule5)
#' }
print.rule <- function(rule) {
  
  #print.rule generates a summary of a rule object
  #it describes if both target outcome and max level can be anything (ie: any hits)
  #it describes if target outcome can be anything (ie: any outcome)
  #as well as the p value
  
  target_outcome = rule$target_outcome
  max_level = rule$max_level
  p_value = rule$p_value
  
  if (target_outcome == "ALL") {
    if (is.null(max_level)) {
      cat("Any hits") #if outcome/max level can be anything
    }
    else {
      cat("Any outcome") #only if outcome can be anything
    }
  }
  else {
    cat(target_outcome, sep = " ")
  }
  
  if (!is.null(max_level)) {
    cat(" with maximum level ")
    cat(max_level)
  }
  
  cat(", p value of ") #add p value at end
  cat(p_value)
  cat("\n")
}