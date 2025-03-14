#' Batch Execution of ROAD Algorithm
#'
#' This function performs the full ROAD workflow, including training, cross-validation, and testing.
#'
#' @param x A numeric matrix of training data (n x p).
#' @param y A numeric vector of training labels (0/1 or 1/2).
#' @param xtest A numeric matrix of test data (n x p).
#' @param ytest A numeric vector of test labels (0/1 or 1/2).
#' @param droad An integer indicating whether to use diagonal ROAD (default: 0).
#' @param sroad An integer indicating whether to use screening ROAD (default: 0).
#' @param para A list of parameters for the algorithm (optional).
#' @return A list containing:
#'   \item{trainError}{Training error rate.}
#'   \item{testError}{Test error rate.}
#'   \item{w}{Optimal weights.}
#'   \item{wPath}{Solution paths for the weights.}
#'   \item{lamList}{List of lambda values.}
#'   \item{num}{Number of non-zero weights.}
#'   \item{lamind}{Index of the optimal lambda.}
#'   \item{cvError}{Cross-validation error rate for the optimal lambda.}
#'   \item{allCvError}{Cross-validation error rates for all lambda values.}
#'   \item{allTestError}{Test error rates for all lambda values.}
#' @export
roadBatch <- function(x, y, xtest, ytest, droad = 0, sroad = 0, para = NULL) {
  
  fit = road(x,y,droad,sroad,para)
  fitCV = roadCV(x,y,fit)
  
  traincl= roadPredict(x, fit, fitCV)
  trainclass=traincl[1][[1]]
  trainclassAll=traincl[2][[1]]
  
  obj<-c()
  obj$trainError = mean(trainclass!=y)
  
  cl = roadPredict(xtest, fit, fitCV)
  class=cl[1][[1]]
  classAll=cl[2][[1]]
  
  obj$testError = mean(class!=ytest)
  obj$w = fitCV$w
  obj$wPath = fit$wPath
  obj$lamList = fit$lamList
  obj$num = fitCV$num
  obj$lamind = fitCV$lamind
  obj$cvError = fitCV$error[fitCV$lamind]
  obj$allCvError =  fitCV$error
  obj$allTestError =  colMeans(classAll!=rep(ytest,fit$para$K))
  
  return(obj)
  
}
