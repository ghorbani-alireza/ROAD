#' Compute Statistics for ROAD Algorithm
#'
#' This function computes necessary statistics (e.g., means, covariance) for the ROAD algorithm.
#'
#' @param x A numeric matrix of predictors (n x p).
#' @param y A numeric vector of labels (0/1 or 1/2).
#' @return A list containing:
#'   \item{x1}{Data for class 1.}
#'   \item{x2}{Data for class 2.}
#'   \item{n1}{Number of samples in class 1.}
#'   \item{n2}{Number of samples in class 2.}
#'   \item{mua}{Mean of the average of class means.}
#'   \item{mud}{Mean difference between classes.}
#'   \item{Sigma}{Pooled covariance matrix.}
#' @export
stat <- function(x, y) {
  if( min(y)==0 ){
    y = y + 1
  }
  x1 = x[y==1,]
  x2 = x[y==2,]
  n1 = dim(x1)[1]
  n2 = dim(x2)[1]
  mu1 = colMeans(x1)
  mu2 = colMeans(x2)
  mua = (mu2+mu1)/2
  mud = (mu2-mu1)/2
  Sigma = (n1*cov(x1) +  n2*cov(x2))/(n1+n2)
  return(list(x1,x2,n1,n2,mua,mud,Sigma))
}
