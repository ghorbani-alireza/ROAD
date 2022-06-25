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