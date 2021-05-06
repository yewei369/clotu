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


suotu<-function(ha){ ## summary of accuracy
  npairs<-length(ha$ana1)
  nms<-NULL
  for(i in 1:npairs) nms<-c(nms,names(ha$ana1)[i])
  max.len0<-max(nchar(nms))
  max.len1<-0
  max.len2<-0
  for(i in 1:length(nms)){
    ch<-nchar(unlist(strsplit(nms[i],"/")))
    if(ch[1]>max.len1) max.len1<-ch[1]
    if(ch[2]>max.len2) max.len2<-ch[2]}
  
  ntotal<-0;acc<-0
  for(i in 1:npairs){
    nm<-names(ha$ana1)[i]
    na<-unlist(strsplit(nm,"/"))
    nm1<-na[1]
    nm2<-na[2]
    
    cat(sprintf("%*s",max.len0+2,"Accuracy"))
    cat(sprintf("%*s",max.len1+2,nm1))
    cat(sprintf("%*s",max.len2+2,nm2))
    cat(sprintf("%*s",7,"Whole"))
    cat("\n")
    cat(sprintf("%*s",max.len0+2,nm))
    cat(sprintf("%*s",max.len1+2,round(mean(ha$ana1[[i]]$accuracy[1,]),2)))
    cat(sprintf("%*s",max.len2+2,round(mean(ha$ana1[[i]]$accuracy[2,]),2)))
    mu<-mean(ha$ana1[[i]]$accuracy[3,])
    cat(sprintf("%*s",7,round(mu,2)))
    cat("\n");cat("\n")
    
    n<-ncol(ha$ana1[[i]][[1]])+ncol(ha$ana1[[i]][[2]])
    ntotal<-ntotal+n;acc=acc+mu*n}
  
  cat("Overall accuracy is ");cat(round(acc/ntotal,2))}