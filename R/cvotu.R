#' Cross validation of OTUs
#' 
#' Implement cross valildation and record the validation error for each iteration and lambda 
#' 
#' @param da, dataset
#' @param ..., mean values by esti.mode=2, available only when esti.mode=2
#' @param lambda, the range of lambda for cross validation, default value of seq(0.001,0.3,by=0.01)
#' @param nfold, number of folds for cross validation, default value of 6
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @return the validation error for each iteration and lambda 
#' 


cvotu<-function(da,...,lambda=seq(0.001,0.3,by=0.01),nfold=6,esti.mode=0){ ## cross validation 
  cver<-rep(0,length(lambda))
  mu<-list(...)[[1]]
  d<-dim(da)[1]
  I<-diag(rep(1,d)) ## entity matrix
  
  trva<-data_split(da,c("train","valid"),c(nfold-1,1))  ## sampling
  ## estimate on train data
  if(esti.mode %in% c(0,1,3)) va_train<-varotu(trva$train,esti.mode=esti.mode)
  else if(esti.mode==2) va_train<-varotu(trva$train,mu,esti.mode=esti.mode)
  
  tryCatch(fit_train<-clime(va_train,lambda,sigma=T),
           error=function(e){fit_train<<-clime::clime(va_train,lambda,sigma=T,perturb=F,linsolver="sim")})
  #fit_train<-clime(va_train,lambda,sigma=T)
  est_omega<-array(0,dim=c(d,d,length(lambda)))
  for(j in 1:length(lambda)) est_omega[,,j]<-fit_train$Omega[[j]]
  ## validate
  
  if(esti.mode %in% c(0,1,3)) va_valid<-varotu(trva$valid,esti.mode=esti.mode)
  else if(esti.mode==2) va_valid<-varotu(trva$valid,mu,esti.mode=esti.mode)
  
  for(j in 1:length(lambda)) cver[j]<-sum(diag(((est_omega[,,j]%*%va_valid)-I)^2)) 
  return (cver)}