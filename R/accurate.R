#' Accurate rate
#' 
#' Calculate accuracy rate for test dataset of a specified target group
#' 
#' @param da, test dataset of a specified target group
#' @param tar, true value for the targets
#' @param sigma, estimated covariance matrix for the given target pair by Clime
#' @param mu1, mean values for the first target group
#' @param mu2, mean values for the second target group
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @param cl.mode, digits for classifying rule, 0 for the one in Kaul's software, 1 for the one in Kaul's paper
#' @return the accuracy for the specific target group under given pair
#' 


accurate<-function(da,tar,sigma,mu1,mu2,esti.mode=0,cl.mode=0){## calculate accuracy rate
  n<-dim(da)[2]
  score<-0
  for(i in 1:n){
    re<-clotu(da[,i],sigma,mu1,mu2,esti.mode,cl.mode)
    if(re==tar) score=score+1}
  return (score/n)}

