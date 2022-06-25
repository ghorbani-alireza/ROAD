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