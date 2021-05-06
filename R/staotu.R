#' Statistics of OTU
#' 
#' Basic statistic values of each OTU
#' 
#' @param da, matrix of OTU 
#' @return a list of mean and variance of those non-zero OTU counts, and index of non-zeros OTUs
#' @examples 
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names
#' staotu(al)


staotu<-function(da){## mean and var for non-zero samples in each OTU
  d<-dim(da)[1];n<-dim(da)[2]
  ind<-sapply(as.data.frame(t(da)),function(x) which(x!=0))
  mu<-rep(NULL,d);va<-rep(NULL,d)
  for(i in 1:d) {mu[i]<-mean(da[i,ind[[i]]])
  va[i]<-var(da[i,ind[[i]]])}
  return (list(mu=mu,var=va,ind=ind))}
