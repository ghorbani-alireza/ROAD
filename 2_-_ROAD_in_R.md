ROAD in R
================

R Code for the high dimensional classification method introduced in the
paper: “A ROAD to Classification in High Dimensional Space (2010)”

In what follows, the code separated into functions:

Stats

``` r
#stat
stat<- function(x, y){
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
```

ROAD Core

``` r
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
```

ROAD

``` r
#road

# This function implements the ROAD
# Algorithm for the following optimization problem:
# minimize 0.5 * w' \Sigma w + \lambda |w|_1 + 1/2 \gamma (Aw-b)^T (Aw-b)
# x: n*p matrix
# y: label for the data, 0/1, or 1/2.

road<-function(x,y,dRoad=0,sRoad=0,para=NULL){
  
  
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
```

Crossvalind

``` r
# Crossvalind

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
```

ROADCV

``` r
#raodCV

# This function implements the ROAD
# Algorithm for the following optimization problem:
# minimize 0.5 * w' \Sigma w + \lambda |w|_1 + 1/2 \gamma (Aw-b)^T (Aw-b)
# x: n*p matrix
# y: label for the data, 0/1, or 1/2.

roadCV<-function(x, y, fit, nfold=5, tracking=0) {

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
```

ROADPredict

``` r
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
```

ROADBatch

``` r
#roadBatch
roadBatch<-function(x, y, xtest, ytest, droad=0, sroad=0, para=NULL){
  
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
```
