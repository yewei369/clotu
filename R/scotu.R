#' Split datasets
#' 
#' Split dataset into training/validation/test sets
#' 
#' @param data,
#' @param labels,
#' @param vec,
#' @return 
#' @examples 
#' 


scotu<-function(da,mu,esti.mode=0){ ## Scale datasets
  d<-dim(da)[1]
  da<-as.matrix(da)
  #re<-matrix(0,nrow=d,ncol=d)
  
  if(esti.mode %in% c(0,1,3)){
    pres<-present(da)
    for(i in 1:d) if(length(pres$ind[[i]])>0) da[i,pres$ind[[i]]]<-da[i,pres$ind[[i]]]-mu[i]  }
  return (da)}