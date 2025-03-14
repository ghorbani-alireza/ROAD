#' Core ROAD Optimization Algorithm
#'
#' This function implements the core optimization algorithm for the ROAD method.
#'
#' @param A A numeric vector representing the mean difference (mud).
#' @param b A numeric scalar (default: 1).
#' @param Sigma A numeric matrix representing the pooled covariance matrix.
#' @param lamList A numeric vector of regularization parameters (lambda values).
#' @param para A list of parameters for the algorithm.
#' @return A list containing:
#'   \item{wPath}{A matrix of solution paths for the weights.}
#'   \item{num}{A vector of the number of non-zero weights for each lambda.}
#' @export
roadCore <- function(A, b, Sigma, lamList, para) {
  
  iterMax = para$iterMax
  epsCon = para$epsCon
  K = para$K
  gamma = para$gamma
  p = para$p
  
  A1List = gamma*A
  cons = A1List^2/gamma + diag(Sigma)
  
  ##   starting the algorithm
  
  wPath = matrix(0,nrow=p,ncol=K)
  w0=matrix(0,nrow=p,ncol=1)
  num = matrix(0,nrow=K,ncol=1)
  pv = 1:p
  
  for (i in 1:K){
    lambda=lamList[i]
    wlast = w0
    
    iterInd = 1
    for (j in 1:p){
      uHat = A1List[j]*b - (Sigma[,j] + gamma*A[j]*A)%*%w0 + (Sigma[j,j]+gamma*A[j]^2)*w0[j]
      w0[j] = sign(uHat)*max(0,(abs(uHat)-lambda))/cons[j]
    }
    error = sum(abs(w0-wlast))
    while(error > epsCon){
      iterInd = 1
      S1 = which(w0!=0)
      
      while( error > epsCon && iterInd < iterMax){
        ###  starting coordinate-wise descent
        iterInd = iterInd + 1
        wlast = w0
        
        for (j in pv[S1]){
          uHat = A1List[j]*b-(Sigma[,j]+gamma*A[j]*A)%*%w0+(Sigma[j,j]+gamma*A[j]^2)*w0[j]
          w0[j] = sign(uHat)*max(0,(abs(uHat)-lambda))/cons[j]
        }
        error = sum(abs(w0-wlast))
      }
      wlast = w0
      for (j in 1:p){
        uHat = A1List[j]*b - (Sigma[,j]+gamma*A[j]*A)%*%w0 + (Sigma[j,j]+gamma*A[j]^2)*w0[j]
        w0[j] = sign(uHat)*max(0,(abs(uHat)-lambda))/cons[j]
      }
      error = sum(abs(w0-wlast))
    }
    if(iterInd == iterMax){
      cat('Unable to converge...at step',i)
    }
    
    wPath[,i] = w0
    num[i] = length(which(w0!=0))
  }
  return(list(wPath,num))
}
