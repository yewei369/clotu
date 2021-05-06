#' Normalizing OTU data
#' 
#' Normalize OTU data with reference OTUs
#' 
#' @param data_tar, dataset split by targets
#' @param thr, threshold for zombie OTU/sample
#' @param method, digits for normalization method, 0 for Kaul's method, 1 for Jun Li's
#' @param ref, index of reference OTUs
#' @param del.otu, default value of T; if T the zombie OTU will be deleted
#' @param del.sam, default value of T; if T the zombie sample will be deleted
#' @return normalized OTU dataset
#' @examples
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234) 
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names
#' ana0<-zomotu(al,0) ## delete otus not present in any sample, no otu names
#' ta<-target_split(da,ana0$otu,"Target")
#' ta_norm<-normalizer(ta,0,0,ana0$ref,FALSE,TRUE)


normalizer<-function(data_tar,thr,method=0,ref,del.otu=TRUE,del.sam=TRUE){

  re<-data_tar
  for(i in 1:length(data_tar)){
    da<-data_tar[[i]]
    d<-dim(da)[1]
    n<-dim(da)[2]
    for(j in 1:d)
      for(k in 1:n)
        if(method==0) re[[i]][j,k]<-ifelse(da[j,k]==0,0,log(da[j,k]/da[ref[1],k]))
    else if(method==1) re[[i]][j,k]<-ifelse(da[j,k]==0,0,log(da[j,k]/mean(da[ref,k])))
    
    
    if(del.otu) {temp<-zomotu(re[[i]],thr)
    if(length(temp$ref)>1) stop(paste("Reference OTU is missing for",i,"-th target!",sep=""))
    re[[i]]<-temp$otu}  ## delete zombie OTUs    
    if(del.sam) {temp<-zomotu(t(re[[i]]),thr)
    if(length(temp$ref)>1) stop(paste("Reference sample is missing for",i,"-th target!",sep=""))
    re[[i]]<-t(temp$otu)} } ## delete zombie samples
  
  return (re)}

