# ROAD Core

# This function implements the ROAD
# Algorithm for the following optimization problem:
# minimize 0.5 * w' \Sigma w + \lambda |w|_1 + 1/2 \gamma (Aw-b)^T (Aw-b)
# x: n*p matrix
# y: label for the data, 0/1, or 1/2.

roadCore<-function(A, b, Sigma, lamList, para){
  
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