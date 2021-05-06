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


simotu.diri40<-function(n,p=693,t=3,ref=c(1,214,490,512,513),seed=F,unif.min=0,unif.max=0.35){ ## simulate OTUs
  ## n, sample size for each target
  ## p, number of OTUs
  ## t, number of targets
  ## this simulation model has extremely high proportion of zeros around 99.5%
  
  meta<-data.frame(SampleID=NULL,Target=NULL)
  otu<-data.frame(OTU.ID=as.factor(c(1:p)))
  if(is.numeric(seed)) set.seed(seed)
  
  mu.sum<-c(40854,43179,47127)
  sd.sum<-c(6987.571,4222.528,6975.728)
  
  #pr<-0
  
  for(i in 1:t){
    da<-pomus$otu[,(36*(i-1)+2):(36*i+1)]
    su<-sapply(as.data.frame(da),function(x) sum(as.numeric(x))) # sum per sample
    den<-as.data.frame(t(sapply(as.data.frame(t(da)),function(x) x/su)))
    #t<<-den
    mu<-as.numeric(sapply(as.data.frame(t(den)),mean))
    var<-as.numeric(sapply(as.data.frame(t(den)),var))
    
    alpha<-rep(0,length=p) ## parameter vector for Dirichlet dist for ith target
    alpha0<-max(mu)*(1-max(mu))/var[which.max(mu)]-1
    alpha<-alpha0*mu
    
    for(j in 1:n){
      d<-runif(1,unif.min,unif.max) ## density of not being structural zeros
      ms<-rbinom(p,1,d)  ## missing structure
      ind<-which(ms==1)
      ind<-union(ind,ref)
      #alpha_sub<-alpha[ind]/sum(alpha[ind])
      
      sc<-round(rnorm(1,mu.sum[i],sd.sum[i]),0) # scale of sample ## mu.sum[i]##
      #prob<-as.vector(rdirichlet(1,alpha_sub))
      prob<-as.vector(rdirichlet(1,alpha))
      prob_sub<-prob[ind]/sum(prob[ind])
      sam<-ceiling(prob_sub*sc) 
      
      sim<-rep(0,p)
      sim[ind]<-sam
      
      na<-paste('sample',nrow(meta)+1,sep="")
      meta<-rbind(meta,data.frame(SampleID=na,Target=paste("target",i,sep="")))
      otu<-cbind(otu,sim)
      names(otu)[ncol(otu)]<-na}    }
  
  meta$SampleID<-as.character(meta$SampleID)
  meta$Target<-as.character(meta$Target)
  return (list(otu=otu,meta=meta))}