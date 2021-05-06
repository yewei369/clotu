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


muotu<-function(da,esti.mode=0){## Mean trva datasets through substracting mean values
  ## da1/da2, train+valid datasets of target groups
  ## mode, 0-Abhishek; 1-JunLi
  d<-dim(da)[1]
  if(esti.mode %in% c(0,1)){
    pres<-present(da)
    mu<-sapply(as.data.frame(t(da)),sum)/(pres$blood*dim(da)[2])
    return (mu)}
  
  else if(esti.mode %in% c(2,3)){
    mu<-array(0,dim=c(2,d,d))
    for(i in 1:d)
      for(j in 1:d){
        ind<-which(da[i,]!=0 & da[j,]!=0)
        mu[,i,j]<-ifelse(length(ind)==0,0,c(sum(da[i,ind])/length(ind)  , sum(da[j,ind])/length(ind)  ) ) } }
  if(esti.mode==2) return (mu)
  else if(esti.mode==3) {
    mu1<-rep(0,d)
    for(i in 1:d) mu1[i]<-mean(mu[1,i,])
    return (mu1)}}  