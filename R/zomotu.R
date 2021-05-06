#' Zombie OTUs
#' 
#' Find zombie OTUs, and preprocess the dataset through deleting zombie OTUs and zombie observations
#' 
#' @param data, matrix of OTU data
#' @param thr, threshold for zombie OTU
#' @param del, when T, delete zombie OTUs
#' @return a list of preprocessed dataset, reference OTUs
#' @examples
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names 
#' zomotu(al,thr=0)



zomotu<-function(data,thr=0,del=TRUE){ ## analyze zombie OTUs, OTU in (otu*sample)
  ##thr, blood threshhold
  ##del, TRUE if zombie OTUs will be deleted
  ##ref, OTUs present in all samples
  
  blood<-present(data)$blood
  nyda<-as.data.frame(data[which(blood>thr),])
  #bacs<-as.character(data$otu[,1])[which(blood>thr)]
  nyblood<-present(nyda)$blood
  ref<-which(nyblood==1)
  
  if(del) return(list(otu=nyda,ref=ref)) else  #,bacs=bacs
    return(list(ref=ref)) }   #,bacs=bacs
