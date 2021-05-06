#' Presency of OTUs
#' 
#' Analyze the presency of OTUs
#' 
#' @param data, matrix of OTU data
#' @return, proportion of non-zero counts for each OTU 
#' @examples
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names 
#' present(al)


present<-function(data){ ## return number of samples each OTU is present in
  d<-dim(data)[1]
  n<-dim(data)[2]
  ind<-list()
  blood<-rep(0,d)
  for(i in 1:d) {
    ind[[i]]<-which(data[i,]!=0)
    blood[i]<-sum(data[i,]!=0)/n}
  return(list(blood=blood,ind=ind))}
