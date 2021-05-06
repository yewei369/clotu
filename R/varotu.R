#' Variance of OTU
#' 
#' Calculate sample covariance matrix for a given target pair
#' 
#' @param da, dataset
#' @param ..., mean values by esti.mode=2, available only when esti.mode=2
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @return covariance matrix of sample OTUs
#' 


varotu<-function(da,...,esti.mode=0){## raw covariance method
  d<-dim(da)[1]
  da<-as.matrix(da)
  re<-matrix(0,nrow=d,ncol=d)
  
  if(esti.mode %in% c(0,1,3)){
    num<-matrix(0,nrow=d,ncol=d) ## number of non-zero samples in both OTUs
    va<-da%*%t(da)
    for(i in 1:d)
      for(j in 1:d){
        num[i,j]<-length(which(da[i,]!=0 & da[j,]!=0))
        re[i,j]<-ifelse(num[i,j]==0,0,va[i,j]/num[i,j])}  }
  else if(esti.mode==2){
    mu<-list(...)[[1]]
    for(i in 1:d)
      for(j in 1:d){
        ind<-which(da[i,]!=0 & da[j,]!=0) ## index of non-zero samples in both OTUs
        va<-t(da[i,ind]-mu[,i,j][1])%*%(da[j,ind]-mu[,i,j][2])
        re[i,j]<-ifelse(length(ind)==0,0,va/length(ind))}    }
  return(re)}