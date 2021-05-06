#' Simulation data with Gaussian distribution
#' 
#' Generate simulation data with Gaussian distribution
#' 
#' @param n, sample size for each target
#' @param p, number of OTUs
#' @param t, number of targets
#' @param nref, number of reference OTUs
#' @param unif.min, min limit of Uniform dist for non-zero probability
#' @param unif.max, max limit of Uniform dist for non-zero probability
#' @param full.mean, full range of mean value of total counts
#' @param full.sd, full range of standard deviation of total counts
#' @param seed, index of the seed, if F (default) no seed is set
#' @return a list of OTU table and meta information
#' @examples 
#' simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)


simotu.gaus<-function(n,p,t,nref=3,unif.min=0,unif.max=0.35,full.mean=20000,full.sd=500,seed=F){ ## simulate OTUs

  meta<-data.frame(SampleID=NULL,Target=NULL)
  otu<-data.frame(OTU.ID=as.factor(c(1:p)))
  if(is.numeric(seed)) set.seed(seed)
  ref<-sample(p,nref)
  
  for(i in 1:t){
    sigma<-matrix(0,nrow=p,ncol=p)
    for(a in 1:p) for(b in 1:p) { #print(0.4+0.1*(i-floor(t/2)))
      sigma[a,b]<-ifelse((0.4+0.1*(i-floor(t/2)))^abs(a-b)<0.01,0,(0.4+0.1*(i-floor(t/2)))^abs(a-b))}
    ref.mean<-sample(full.mean,p,replace=T)
    ref.sd<-sample(full.sd,p,replace=T)
    
    for(j in 1:n){
      d<-runif(1,unif.min,unif.max)
      ms<-rbinom(p,1,d)
      ind<-which(ms==1)
      ind<-union(ind,ref)
      sub_sigma<-sigma[ind,ind]
      if(length(ind)>0){
        if(length(ind)==1) sam<-rnorm(1,mean=rep(0,length(ind)),sd=sqrt(sub_sigma)) else
          sam<-mvtnorm::rmvnorm(1,mean=rep(0,length(ind)),sigma=sub_sigma)
        
        sam<-round(sam*ref.sd[ind]+ref.mean[ind],0)
        #if(min(sam)<0) sam<-sam-min(sam)+1
        sam[which(sam==0)]<-1
        
        sim<-rep(0,p)
        sim[ind]<-sam} else
          sim<-rep(0,p)
      
      na<-paste('sample',nrow(meta)+1,sep="")
      meta<-rbind(meta,data.frame(SampleID=na,Target=paste("target",i,sep="")))
      otu<-cbind(otu,sim)
      names(otu)[ncol(otu)]<-na}    }
  for(j in 2:(t*n+1))
    if(length(which(otu[j,-1]<0))>0) 
      otu[j,-1][which(otu[j,-1]!=0)]<-otu[j,-1][which(otu[j,-1]!=0)]-min(otu[j,-1][which(otu[j,-1]!=0)])+1
  
  meta$SampleID<-as.character(meta$SampleID)
  meta$Target<-as.character(meta$Target)
  return (list(otu=otu,meta=meta))}