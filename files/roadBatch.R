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
