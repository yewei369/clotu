#' Split dataset
#' 
#' Split dataset into training/validation/test sets
#' 
#' @param data, matrix of data
#' @param labels, character vector of labels for splitted datasets
#' @param vec, proportion of the splitted datasets
#' @return list of splitted datasets
#' @examples 
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234)  
#' al<-data_extract(da,Target %in% c("target1","target2","target3")) # no otu names
#' data_split(al,c("train","valid","test"),c(0.6,0.2,0.2))


data_split<-function(data,labels,vec){  
  
  re<-list()
  n<-dim(data)[2]
  nm<-colnames(data)
  ind_all<-1:n
  
  for(i in 1:length(labels)){
    ind<-sample(n,round(n*vec[i]/sum(vec)))
    re[[labels[i]]]<-as.data.frame(data[,ind])
    names(re[[labels[i]]])<-nm[ind]}
  
  return (re)}
