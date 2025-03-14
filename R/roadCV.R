#' Cross-Validation for ROAD Algorithm
#'
#' This function performs k-fold cross-validation for the ROAD algorithm.
#'
#' @param x A numeric matrix of predictors (n x p).
#' @param y A numeric vector of labels (0/1 or 1/2).
#' @param fit A fitted ROAD model object.
#' @param nfold An integer specifying the number of folds for cross-validation (default: 5).
#' @param tracking An integer indicating whether to display progress (default: 0).
#' @return A list containing:
#'   \item{error}{Cross-validation error rates.}
#'   \item{lamind}{Index of the optimal lambda.}
#'   \item{w}{Optimal weights.}
#'   \item{num}{Number of non-zero weights.}
#' @export
roadCV <- function(x, y, fit, nfold = 5, tracking = 0) {
  
  ###    K: number of lambdas on the log-scale
  ###    epsilon: lamMin=lamMax*epsilon
  
  para = fit$para
  
  if(!is.null(para$nfold)){
    nfold=para$nfold
  }
  
  K = para$K
  
  dRoad = para$dRoad
  if(para$sRoad!=0){
    x = x[,para$sInd]
  }
  
  b = 1
  n=dim(x)[1]
  p=dim(x)[2]
  
  if(min(y)==0){
    y = y + 1
  }
  
  x1 = x[y==1,]
  x2 = x[y==2,]
  n1 = dim(x1)[1]
  n2 = n-n1;
  
  indices1  = crossvalind(n1, nfold)
  indices2  = crossvalind(n2, nfold)
  
  errorMat = matrix(0,nfold, K)
  
  SumDist = matrix(0,1,K)
  SumDist2= SumDist
  
  for (i in 1:nfold){
    if(tracking==1){
      cat('CV ', i, ' of ', nfold)
    }
    
    test1 = (indices1==i)
    test2 = (indices2==i)
    train1 = !test1
    train2 = !test2
    x1train = x1[train1,]
    x2train = x2[train2,]
    x1test = x1[test1,]
    x2test = x2[test2,]
    n1test = dim(x1test)[1]
    n2test = dim(x2test)[1]
    n1train = n1-n1test
    n2train = n2-n2test
    mu1 = colMeans(x1train)
    mu2 = colMeans(x2train)
    mua = (mu2+mu1)/2
    mud = (mu2-mu1)/2
    Sigma = (n1train*cov(x1train) + n2train*cov(x2train))/(n1train+n2train)
    if(dRoad!=0){Sigma = diag(diag(Sigma))}
    
    A = mud
    
    wn=roadCore(A, b, Sigma, fit$lamList, para)
    wPath=wn[1][[1]]
    num=wn[2][[1]]
    
    x1T = (x1test-matrix(rep(mua, n1test),n1test,byrow=T))
    x2T = (x2test-matrix(rep(mua, n2test),n2test,byrow=T))
    
    x1Tnorm = sqrt(colSums(x1T^2))
    x2Tnorm = sqrt(colSums(x2T^2))
    
    x1TnormMat = matrix(rep(x1Tnorm, n1test),nc=n1test)
    x2TnormMat = matrix(rep(x2Tnorm, n2test),nc=n2test)
    
    dist1 = (x1test-matrix(rep(mua, n1test),n1test,byrow=T))%*%wPath
    dist2 = -(x2test-matrix(rep(mua, n2test),n2test,byrow=T))%*%wPath
    test1 = dist1 >0
    test2 = dist2 >0
    
    wPathnorm = sqrt(colSums(wPath^2))
    wPathnormmat1 = matrix(rep(wPathnorm, n1test),n1test,byrow=T)
    wPathnormmat2 = matrix(rep(wPathnorm, n2test),n2test,byrow=T)
    dist1 = dist1/wPathnormmat1
    dist2 = dist2/wPathnormmat2
    dist3 = dist1^2
    dist4 = dist2^2
    
    distAve = colMeans(rbind(dist1, dist2))
    distAve2 = colMeans(rbind(dist3, dist4))
    
    SumDist = SumDist + distAve
    SumDist2 = SumDist2 + distAve2
    
    errorMat[i,] = errorMat[i,]+ colSums(rbind(test1,test2))
  }
  
  fitCV = c()
  fitCV$error = colSums(errorMat)/n
  fitCV$lamind = max(which(fitCV$error==min(fitCV$error)))
  fitCV$w = fit$wPath[, fitCV$lamind]
  fitCV$num = sum(fitCV$w!=0)
  
  return(fitCV)
}
