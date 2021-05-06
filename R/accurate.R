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


accurate<-function(da,tar,sigma,mu1,mu2,esti.mode=0,cl.mode=0){## calculate accuracy rate
  n<-dim(da)[2]
  score<-0
  for(i in 1:n){
    re<-clotu(da[,i],sigma,mu1,mu2,esti.mode,cl.mode)
    if(re==tar) score=score+1}
  return (score/n)}

