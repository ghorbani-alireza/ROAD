#' Generate Cross-Validation Indices
#'
#' This function generates indices for k-fold cross-validation.
#'
#' @param N An integer specifying the total number of samples.
#' @param kfold An integer specifying the number of folds.
#' @return A numeric vector of indices for cross-validation.
#' @export
crossvalind <- function(N, kfold) {
  n1 = N%%kfold
  a = sample(1:kfold,n1)
  n2 = N%/%kfold
  for(i in 1:n2){
    a = c(a,sample(1:kfold,kfold))
    a = sample(a,length(a))
  }
  return(a)
} 
