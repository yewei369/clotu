#' Major rules of classification
#' 
#' Implement major rules to classify an arbitrary observation against all targets
#' 
#' @param tr, trained object by trotu function
#' @param shrink.mode, digits for shrinking method, 0 for Kaul's method, 1 for Jun Li's
#' @param esti.mode, digits for covariance estimating method, 0 for Kaul's method, 1/2/3 for Jun Li's 1st/2nd/3rd method
#' @param cl.mode, digits for classifying rule, 0 for the one in Kaul's software, 1 for the one in Kaul's paper
#' @param thr, threshold for zombie OTU
#' @param nvar, number of OTUs after shrinking
#' @return Table of winning targets and its winning times

#da<-simotu.gaus(50,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234) 
#ha<-trotu(da,Target %in% c("target1","target2","target3"),
#           thr=0,target="Target",pairs=c("target1","target2"),del.otu=FALSE,del.sam=TRUE,nvar=75,
#           lambda=seq(0.001,0.3,by=0.01),nsim=3,seed=FALSE,nfold=5, nsampling=1, 
#           test.per=0.2,norm.mode=1,shrink.mode=1,esti.mode=2,cl.mode=0)
#majorcl(ha,1,2)


majorcl<-function(tr,shrink.mode,esti.mode,cl.mode=0,thr=0,nvar=75){
  tic<-Sys.time()
  pre<-NULL
  tru<-NULL
  
  
  da<-tr$ta_norm
  ntar<-length(da)
  nmtar<-names(da)
  npair<-length(tr$ana1)
  nmpairs<-names(tr$ana1)
  nsim<-dim(tr$ana1[[1]]$sigma)[3]
  
  board<-array(0,dim=c(npair,3,dim(tr$ana0$otu)[2]))
  rank<-list()
  pointer<-0
  for(i in 1:ntar){
    truetar<-nmtar[i]  ## true taret value
    n<-dim(da[[i]])[2] ## sample size in this target
    for(j in 1:n){
      tru<-c(tru,truetar)
      ob<-da[[i]][,j] ## the observation to classify
      pointer<-pointer+1
      
      for(m in 1:npair){
        nmpair<-nmpairs[m]
        allcha<-unlist(strsplit(nmpair,""))
        speind<-which(allcha=="/")
        nmtar1<-substr(nmpair,1,speind-1)
        nmtar2<-substr(nmpair,speind+1,nchar(nmpair))
        
        shda<-shrink_var(da[[nmtar1]],da[[nmtar2]],thr,nvar,c(nmtar1,nmtar2),shrink.mode)
        shob<-ob[shda$ind]
        #print(shda[[1]][1:5,1:5])
        
        if(esti.mode %in% c(0,1)){
          mu1<-muotu(shda[[1]],esti.mode);mu2<-muotu(shda[[2]],esti.mode)## mean values of trva from 2 target datasets
          pa1<-mu1;pa2<-mu2 }## parameters for classification
        else if(esti.mode %in% c(2,3)){
          pa1<-muotu(shda[[1]],esti.mode=3)
          pa2<-muotu(shda[[2]],esti.mode=3)}
        
        count1<-0;count2<-0
        for(n in 1:nsim){
          re<-clotu(shob,tr$ana1[[m]]$sigma[,,n],pa1,pa2,esti.mode,cl.mode)
          if(re==1) count1<-count1+1 else
            count2<-count2+1}
        if(count1>count2) {
          board[m,1,pointer]<-nmtar1
          board[m,2,pointer]<-nmtar2
          board[m,3,pointer]<-count1/nsim} else {
            board[m,1,pointer]<-nmtar2
            board[m,2,pointer]<-nmtar1
            board[m,3,pointer]<-count2/nsim}
      }
      
      ## implement majority rule
      temp<-data.frame()  ## 1st col:winning target; 2nd col:winning times; 3rd col: avg success rate 
      for(m in 1:npair){
        if(m==1) {
          temp<-rbind(temp,data.frame(board[m,1,pointer],1.00,as.numeric(board[m,3,pointer])))
        } else {
          ind<-which(temp[,1]==board[m,1,pointer])
          if(length(ind)>0) {
            temp[ind,2]<-temp[ind,2]+1
            temp[ind,3]<-(temp[ind,3]*(temp[ind,2]-1)+as.numeric(board[m,3,pointer]))/temp[ind,2]
          } else {
            temp<-rbind(temp,data.frame(board[m,1,pointer],1,as.numeric(board[m,3,pointer])))
          }
        }
      }
      
      pre<-c(pre,as.character(temp[order(temp[,2],temp[,3],decreasing=T),][1,1]))
      rank<-append(rank,list(temp))
      
    }  }
  ct<-table(pre,tru)
  print(ct)
  print(paste("Classification accuracy with Majority rule is:",sum(diag(ct))/sum(ct) ))
  
  toc<-Sys.time()
  cat("DONE! | Total time cost: ",round(difftime(toc,tic,units="secs")[[1]]/60,2)," min",sep="")
  return (list(pre=pre,tru=tru,board=board,rank=rank))}
