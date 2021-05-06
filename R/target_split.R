#' Split into target groups
#' 
#' Split extracted data into different target groups
#' 
#' @param data_ori, original dataset consisting dataframes of "otu" and "meta"
#' @param data_pro, extracted matrix of specified condition
#' @param target, target variable
#' @return a list of dataframes for each target
#' @examples 
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names
#' ana0<-zomotu(al,0) ## delete otus not present in any sample, no otu names
#' ta<-target_split(da,ana0$otu,"Target")



target_split<-function(data_ori,data_pro,target){
  ## split processed dataset into different target datasets
  
  targets<-unique( eval(parse(text=paste(data_ori,"$meta$",target,sep="")))) #deparse(substitute(
  target_nr<-length(targets ) ## number of unique targets
  sample_nr<-dim(data_pro)[2] ## number of samples
  samples<-colnames(data_pro) ## SampleID
  re<-list()
  for(i in 1:sample_nr){
    tar<-eval(parse(text=paste("dplyr::filter(data_ori$meta,SampleID==samples[i])$",target,sep="")))  ## target for this sample
    tar<-unlist(strsplit(tar," "))
    if(length(tar)>1) tar<-paste(tar,collapse="_")
    temp<-as.data.frame(eval(parse(text=paste("data_pro$",samples[i],sep=""))))
    names(temp)<-samples[i]
    if(!(tar %in% names(re))) 
    {eval(parse(text=paste("re$",tar,"<-temp",sep="")))}  else
      eval(parse(text=paste("re$",tar,"<-cbind(","re$",tar,",temp)",sep="")))
  }
  return (otu=re)}
