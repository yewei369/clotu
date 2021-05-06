#' Extract data
#' 
#' Extract OTU data for specified condition
#' 
#' @param data, a list consisting of 2 data frames: otu, each column is an observation and each row is a serie of OTU counts for all observations; meta, each row contains all information for that observation 
#' @param ..., filter condition for special subset of the given data
#' @return a matrix of OTUs meeting the specified condition
#' @examples
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' ha<-data_extract(da,Target %in% c("target1","target2","target3")) ## return a matrix with 693 OTUs and 18 observations


data_extract<-function(data,...){
  return(sapply(data$otu[,dplyr::filter(data$meta,...)$SampleID],as.numeric))}