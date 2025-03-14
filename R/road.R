#' ROAD Algorithm for High-Dimensional Classification
#'
#' This function implements the ROAD algorithm for high-dimensional classification.
#'
#' @param x A numeric matrix of predictors (n x p).
#' @param y A numeric vector of labels (0/1 or 1/2).
#' @param dRoad Indicator for diagonal ROAD (default: 0).
#' @param sRoad Indicator for screening ROAD (default: 0).
#' @param para A list of parameters (optional).
#' @return A list containing the fitted model, error rates, and other outputs.
#' @export
road <- function(x, y, dRoad = 0, sRoad = 0, para = NULL) {
  
  
  ######    dRoad: indicator of doing diagonal-road or not.
  ######    sRoad: indicator of doing screening-road or not. sRoad=1: first
  ######    version, sRoad=2: second version.
  ######    K: number of lambdas on the log-scale
  ######    epsilon: lamMin=lamMax*epsilon
  
  if (is.null(para)) {
    para$alpha=0
    para$iterMax=100
    para$epsCon=1e-2
    para$epsilon=1e-3
    para$K=100
    para$gamma=10
  }
  
  
  iterMax = para$iterMax
  epsCon = para$epsCon
  epsilon = para$epsilon
  K = para$K
  gamma = para$gamma
  alpha = para$alpha
  para$dRoad = dRoad
  para$sRoad = sRoad
  
  
  b = 1
  n = dim(x)[1]
  p = dim(x)[2]
  
  para$p = p
  
  xx = stat(x,y)
  x1 = xx[1][[1]]
  x2 = xx[2][[1]]
  n1 = xx[3][[1]]
  n2 = xx[4][[1]]
  mua = xx[5][[1]]
  mud = xx[6][[1]]
  Sigma = xx[7][[1]]
  
  if(dRoad!=0){
    Sigma = diag(diag(Sigma))
  }
  if(sRoad!=0){
    
    Te = mud/sqrt(diag(Sigma))
    per = sample(n)
    yper = y[per]
    
    statper = stat(x,yper)
    mua = statper[5][[1]]
    mudPer = statper[6][[1]]
    SigmaPer = statper[7][[1]]
    
    Tper = mudPer/sqrt(diag(SigmaPer))
    sInd = which(abs(Te)>quantile(abs(Tper),1-alpha))
    
    if(sRoad==2){
      n0 = length(sInd)
      if(2*n0 >= p){
        sInd = (1:p)
      }else {
        index2 = sInd
        index3 = sInd
        tmpMat = abs(Sigma[index2,setdiff(1:p,index2)])
        index4 = setdiff(1:p,index2)
        
        tmpvec = tmpMat
        dim(tmpvec) = c(n0*(p-n0),1)
        
        for (i in 1:n0){
          tmpvecmax=max(tmpvec)
          ind=which.max(tmpvec)
          aa1 = ceiling(ind/n0)
          aa2 = ind%%n0
          
          if(aa2==0){
            aa2=n0
          }
          
          dim1 = index4[aa1]
          index3[i] = dim1
          tmpvec[seq((aa2-1)+1,(n0*(p-n0)),n0)]=0
          tmpvec[n0*(aa1-1)+(1:n0)]=0
        }
        sInd = c(index2, index3)
      }
    }
    
    para$sInd = sInd
    para$p = length(sInd)
    x = x[,sInd]
    xx = stat(x,y)
    x1 = xx[1][[1]]
    x2 = xx[2][[1]]
    n1 = xx[3][[1]]
    n2 = xx[4][[1]]
    mua = xx[5][[1]]
    mud = xx[6][[1]]
    Sigma = xx[7][[1]]
    
  }
  
  lamMax = max(abs(gamma*mud))*(1-epsilon)
  lamMin = lamMax*epsilon
  lamList =logspace(log10(lamMax),log10(lamMin),n=K)
  
  A = mud
  wn = roadCore(A, b, Sigma, lamList, para)
  wPath = wn[1][[1]]
  num = wn[2][[1]]
  
  test1 = (x1-rep(mua, n1))%*%wPath >0 
  test2 = (x2-rep(mua, n2))%*%wPath <0 
  error = colMeans(rbind(test1,test2))
  
  fit = c()
  fit$lamList = lamList
  fit$wPath = wPath
  fit$num = num
  fit$error = error
  fit$para = para
  fit$mua = mua
  
  return(fit)
}
