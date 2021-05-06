#' Classify OTU dataset 
#' 
#' For each target variable pair, divide the given OTU dataset into "trva" (including training and validation set) and test sets, then find an optimal hyperparameter lambda through cross validation, estimate the covariance matrix from "trva" set and obtain the accuracy by classifying the test set.     
#' 
#' @param data, a list consisting of 2 data frames: otu, each column is an observation and each row is a serie of OTU counts for all observations; meta, each row contains all information for that observation 
#' @param ..., filter condition for special subset of the given data
#' @param thr, threshold for zombie OTU
#' @param target, name for target variable  
#' @param pairs, character vector of specified target pairs, should be of length 2; otherwise if F (default value), it will classify all possible target pairs
#' @param del.otu, default value of T; if T the zombie OTU will be deleted
#' @param del.sam, default value of T; if T the zombie sample will be deleted
#' @param nvar, number of OTUs after shrinking
#' @param lambda, the range of lambda for cross validation, default value of seq(0.001,0.3,by=0.01)
#' @param nsim, number of simulations, default value of 25
#' @param seed, index of the seed, if F (default) no seed is set
#' @param nfold, number of folds for cross validation, default value of 6
#' @param nsampling, number of iterations for cross validation, default value of 20
#' @param test.per, percentage of test dataset, default value of 0.2
#' @param norm.mode, digits for normalization method, 0 for Kaul's method, 1 for Jun Li's
#' @param shrink.mode, digits for shrinking method, 0 for Kaul's method, 1 for Jun Li's
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @param cl.mode, digits for classifying rule, 0 for the one in Kaul's software, 1 for the one in Kaul's paper
#' @return  a list consisting of following elements:
#' ana0, analysis result after deleting zombie OTUs. A list of generated OTU tabel, reference OTUs and the name of OTUs
#' ta, split the OTU table in "ana0" into different targets
#' ta_norm, normalized OTU tables
#' targets, vector of target names
#' vs, pairs of targets
#' ana1, a list of shrinked OTU table for each target, estimated covariance matrix "Sigma", estimated precision matrix "Omega", cross validation error, and accuracies
#' @references 
#' A. Kaul, O. Davidov and S. D. Peddada, "Structural zeroes in high-dimensional data with applications to microbiome studies", Biostatistics, vol. 18, no. 3, p. 422-433, 2017. 
#' Jun Li, "Classification of microbiome data with structural zeroes and small samples", master thesis at Link\"{o}ping University, 2021
#' @examples
#' da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234) 
#' ha<-trotu(da,Target %in% c("target1","target2","target3"),
#'           thr=0,target="Target",pairs=c("target1","target2","target3"),del.otu=F,del.sam=T,nvar=75,
#'           lambda=seq(0.001,0.3,by=0.01),nsim=3,seed=F,nfold=5, nsampling=1, 
#'           test.per=0.2,norm.mode=1,shrink.mode=1,esti.mode=2,cl.mode=0)
#'



trotu<-function(data,...,thr=0,target,pairs=F,
                del.otu=T,del.sam=T,nvar=F,
                lambda=seq(0.001,0.3,by=0.01),nsim=25,
                seed=F,nfold=6, nsampling=20, test.per=0.2,
                norm.mode=0,shrink.mode=0,esti.mode=0,cl.mode=0){

  tic<-Sys.time()
  
  
  al<-data_extract(data,...) # no otu names
  ana0<-zomotu(al,thr) ## delete otus not present in any sample, no otu names
  ta<-target_split(data,ana0$otu,target) ## split into OTUs for different target groups
  
  if(length(ana0$ref)<1) stop("Reference OTU is missing!")
  ta_norm<-normalizer(ta,thr,norm.mode,ana0$ref,del.otu,del.sam) 
  
  targets<-names(ta_norm) ## names of targets/rivers
  tar_nr<-length(ta_norm) ## number of targets
  
  
  if(is.numeric(seed)) set.seed(seed) ## set seed
  pair_pool<-combn(tar_nr,2)  ## pool of pairs, 2*npair
  npair<-ncol(pair_pool)
  vs<-NULL  ## target pairs
  if(is.character(pairs) & length(pairs)==2){
    vs<-matrix(match.arg(pairs,targets,several.ok=T),ncol=2)} else {
      for(i in 1:npair) vs<-rbind(vs,targets[pair_pool[,i]])}
  
  npair<-nrow(vs)  ## actual number of pairs
  ana1<-list()  ## analysis of target pairs
  
  
  for(p in 1:npair){ ## loop each pair of target
    vs_pair<-vs[p,] ## pair
    ## reduce number of variables 
    pair_name<-paste(vs_pair[1],"/",vs_pair[2],sep="")
    ana1[[pair_name]]<-shrink_var(ta_norm[[vs_pair[1]]],ta_norm[[vs_pair[2]]],thr,nvar,vs_pair,shrink.mode)
    d<-dim(ana1[[pair_name]][[1]])[1]  ## OTU number
    
    
    ## validate, train and test model on given data
    
    omega<-array(0,dim=c(d,d,nsim))
    sigma<-array(0,dim=c(d,d,nsim))
    cver<-array(0,dim=c(nsampling,length(lambda),nsim))
    accuracy<-array(0,dim=c(3,nsim))
    
    
    for(s in 1:nsim){ ## generating different training/valid/test datasets
      
      
      tictic<-Sys.time()
      
      da1<-data_split(ana1[[pair_name]][[1]],c("trva","test"),c(1-test.per,test.per))  ## sampling
      da2<-data_split(ana1[[pair_name]][[2]],c("trva","test"),c(1-test.per,test.per))
      da<-cbind(da1$trva,da2$trva)
      if(esti.mode %in% c(0,1)){
        mu1<-muotu(da1$trva,esti.mode);mu2<-muotu(da2$trva,esti.mode)## mean values of trva from 2 target datasets
        pa1<-mu1;pa2<-mu2 }## parameters for classification
      else if(esti.mode %in% c(2,3)){
        pa1<-muotu(da1$trva,esti.mode=3)
        pa2<-muotu(da2$trva,esti.mode=3)}
      
      for(i in 1:nsampling){
        ## cross validation based on target pairs
        ## find optimal Lambda, 
        
        if(esti.mode==0) { ## Abhishek: mean value based on each target
          sc1<-scotu(da1$trva,mu1,esti.mode);sc2<-scotu(da2$trva,mu2,esti.mode)
          cver[i,,s]<-cvotu(cbind(sc1,sc2),lambda,nfold,d,I,esti.mode) }
        else if(esti.mode==1) {## mean value based on whole trva dataset
          mu<-muotu(da,esti.mode) 
          sc<-scotu(da,mu,esti.mode)
          cver[i,,s]<-cvotu(sc,lambda,nfold,d,I,esti.mode) }
        else if(esti.mode==2) {## customized mean, bilateral variance
          mu<-muotu(da,esti.mode) 
          cver[i,,s]<-cvotu(da,mu,lambda,nfold,d,I,esti.mode) }
        else if(esti.mode==3){## blanket mean values
          mu<-muotu(da,esti.mode) 
          sc<-scotu(da,mu,esti.mode)
          cver[i,,s]<-cvotu(sc,lambda,nfold,d,I,esti.mode) }
      }
      

      mu_cver<-rep(0,length(lambda))
      for(i in 1:length(lambda)) mu_cver[i]<-mean(cver[,i,s])
      best_lambda<-which.min(mu_cver)
      
      ## precision estimates, estimate Sigma+Omega matrix
      if(esti.mode %in% c(0,1,3)) va<-varotu(da,esti.mode) 
      else if(esti.mode==2) va<-varotu(da,mu,esti.mode)
      tryCatch(fit<-clime(va,lambda,sigma=T),
               error=function(e){fit<<-clime(va,lambda,sigma=T,perturb=F,linsolver="sim")})
      omega[,,s]<-fit$Omega[[best_lambda]]
      sigma[,,s]<-solve(omega[,,s])
      
      ## classify
      ### accuracy of classifying the first/second/whole target in pair
      n1<-dim(da1$test)[2];n2<-dim(da2$test)[2]
      ac1<-accurate(da1$test,1,sigma[,,s],pa1,pa2,esti.mode,cl.mode) 
      ac2<-accurate(da2$test,2,sigma[,,s],pa1,pa2,esti.mode,cl.mode)
      ac<-ac1*n1/(n1+n2)+ac2*n2/(n1+n2)
      accuracy[,s]<-c(ac1,ac2,ac)
      
      toctoc<-Sys.time()
      cat(paste("Target pairs: ",pair_name," | Simulation: ",s,"/",nsim," | Time cost: ",round(difftime(toctoc,tictic,units="secs")[[1]]/60,2)," min",sep="")) ;   cat("\n")  }
    
    
    
    ana1[[pair_name]]$omega<-omega
    ana1[[pair_name]]$sigma<-sigma
    ana1[[pair_name]]$cver<-cver
    ana1[[pair_name]]$accuracy<-accuracy}
  
  
  ## anao, prelilminary analysis on raw data
  ## ta_norm, normalized targeted datasets
  ## sim, different splitings of normalized targeted datasets 
  toc<-Sys.time()
  cat("DONE! | Total time cost: ",round(difftime(toc,tic,units="secs")[[1]]/60,2)," min",sep="")
  
  return (list(ana0=ana0,ta=ta,ta_norm=ta_norm,targets=targets,
               vs=vs,ana1=ana1))}

