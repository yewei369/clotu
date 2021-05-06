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


clotu<-function(da,sigma,mu1,mu2,esti.mode=0,cl.mode=0){  ## classification
  ##cl.mode: 0-formular in code; 1-formular in paper
  ind<-which(da!=0)
  va<-sigma[ind,ind] ## variance for non-zero OTUs
  tom<<-va
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
      return(sample(2,1))})
}  
## clotu(da1$test[,1],)