#' Classification of one OTU observation
#' 
#' Classify the given observation against a given target pair
#' 
#' @param da, an observation
#' @param sigma, estimated covariance matrix for the given target pair by Clime
#' @param mu1, mean values for the first target group
#' @param mu2, mean values for the second target group
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @param cl.mode, digits for classifying rule, 0 for the one in Kaul's software, 1 for the one in Kaul's paper
#' @return predicted value of the given data
#' 


clotu<-function(da,sigma,mu1,mu2,esti.mode=0,cl.mode=0){  ## classification
  ##cl.mode: 0-formular in code; 1-formular in paper
  ind<-which(da!=0)
  
  if(length(ind)>0){
  va<-sigma[ind,ind] ## variance for non-zero OTUs
  prec<-solve(va) ## precision
  
  sc1<-da[ind]-mu1[ind];sc2<-da[ind]-mu2[ind]
  
  if(cl.mode==0) {cl1<--t(sc1)%*%prec%*%sc1+log(det(prec))
  cl2<--t(sc2)%*%prec%*%sc2+log(det(prec))}
  else if(cl.mode==1) {cl1<-t(da[ind])%*%prec%*%sc1-0.5*t(sc1)%*%prec%*%sc1
  cl2<-t(da[ind])%*%prec%*%sc2-0.5*t(sc1)%*%prec%*%sc2}
  #print(cl1);print(cl2)
  tryCatch({if(cl1>cl2) return(1) else
    return (2)},
    error=function(e) {#message(paste("Comparing between ",cl1," and ",cl2," is not valid!",sep=""))
      return(sample(2,1))})  } else {
        return(sample(2,1))}
}  
