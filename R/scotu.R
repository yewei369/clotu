#' Scaling OTUs
#' 
#' Scale dataset by substructing mean values
#' 
#' @param da, dataset
#' @param mu, mean values of OTUs
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @return scaled counts for all observations
#' 


scotu<-function(da,mu,esti.mode=0){ ## Scale datasets
  d<-dim(da)[1]
  da<-as.matrix(da)
  #re<-matrix(0,nrow=d,ncol=d)
  
  if(esti.mode %in% c(0,1,3)){
    pres<-present(da)
    for(i in 1:d) if(length(pres$ind[[i]])>0) da[i,pres$ind[[i]]]<-da[i,pres$ind[[i]]]-mu[i]  }
  return (da)}