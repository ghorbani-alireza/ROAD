# roadPredict

# This function implements the ROAD
# Algorithm for the following optimization problem:
# minimize 0.5 * w' \Sigma w + \lambda |w|_1 + 1/2 \gamma (Aw-b)^T (Aw-b)
# x: n*p matrix
# y: label for the data, 0/1, or 1/2.

roadPredict<-function(xtest, fit, fitCV){
  
  ###    K: number of lambdas on the log-scale
  ###    epsilon: lamMin=lamMax*epsilon
  
  w = fitCV$w
  mua = fit$mua
  wPath = fit$wPath
  
  ntest = dim(xtest)[1]
  if (fit$para$sRoad!=0){
    xtest = xtest[,fit$para$sInd]
  }
  class = (xtest-matrix(rep(mua, ntest),ntest,byrow=T))%*%w>0
  classAll = (xtest-matrix(rep(mua, ntest),ntest,byrow=T))%*%wPath>0
  
  return(list(class,classAll))
}
