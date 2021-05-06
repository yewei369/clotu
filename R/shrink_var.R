#' Split datasets
#' 
#' Split dataset into training/validation/test sets
#' 
#' @param da1, normalized dataset for target 1
#' @param da2, normalized dataset for target 2
#' @param thr, threshold for zombie OTU
#' @param nvar, number of OTUs after shrinking
#' @param pair, character vector of specified target pairs, should be of length 2
#' @param method, digits for shrinking method, 0 for Kaul's method, 1 for Jun Li's
#' @return a list of shrinked OTU tables for each target
#' @examples
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names
#' ana0<-zomotu(al,0) ## delete otus not present in any sample, no otu names
#' ta<-target_split(da,ana0$otu,"Target")
#' ta_norm<-normalizer(ta,0,1,ana0$ref,FALSE,TRUE)
#' shrink_var(ta_norm[[1]],ta_norm[[2]],0,75,method=1)


shrink_var<-function(da1,da2,thr,nvar,pair=c(1,2),method=0){## shrink number of variables/OTUs
  ## method: 0 for Abhishek method, ie 
  ##         1 for JunLi method, ie 
  sta1<-staotu(da1);sta2<-staotu(da2)
  n1<-length(sta1$ind);n2<-length(sta2$ind)
  if(method==0){
    df<-((sta1$var/n1)+(sta2$var/n2))^2/(((sta1$var/n1)/(n1-1))+((sta2$var/n2)/(n2-1)))} 
  else if(method==1) df<-n1+n2-2
  
  sta<-(sta1$mu-sta2$mu)/sqrt((((n1-1)*sta1$var+(n2-1)*sta2$var)/(n1+n2-2)) *((1/n1)+(1/n2)) )  
  t=2*pt(abs(sta), df, lower.tail=FALSE)
  signif_ind<-order(t)[1:nvar]
  ##### added after simulation test
  re<-list(t(zomotu(t(da1[signif_ind,]),thr)$otu),t(zomotu(t(da2[signif_ind,]),thr)$otu ))
  #re<-list(da1[signif_ind,],da2[signif_ind,])
  names(re)<-pair
  
  return (re)}