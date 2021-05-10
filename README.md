clotu is a public repository developped from an original master thesis "Classification of microbiome data with structural zeroes and small samples"

To reproduce the results in section 5:
1. Classification of Mussel data (run the following code with corresponding digits of norm.mode/shrink.mode/esti.mode/cl.mode, "pomus" is the mussel data which due to agreement with data owner is not public available yet)

   tr<-trotu(pomus,tissue %in% c("gut","foot","gonads","muscles","mantle","gills"),thr=0,target="river",
         pairs=c("C","P","S"),del.otu=F,del.sam=T,nvar=75,lambda=seq(0.001,0.3,by=0.01),nsim=50,seed=F,
         nfold=5,nsampling=20,test.per=0.2,norm.mode=0,shrink.mode=0,esti.mode=2,cl.mode=0)
   majorcl(tr,shrink.mode=0,esti.mode=2)
         
2. Classification of simulation data by Gaussian distribution
   da<-simotu.gaus(num,700,3,nref=5,full.mean=10000,unif.min=0,unif.max=0.4,seed=1234) ##num is sample size for each target
   tr<-trotu(da,Target %in% c("target1","target2","target3"),thr=0,target="Target",
         pairs=c("target1","target2","target3"),del.otu=F,del.sam=T,nvar=75,
         lambda=seq(0.001,0.3,by=0.01),nsim=50,seed=F,nfold=5, nsampling=20, 
         test.per=0.2,norm.mode=1,shrink.mode=1,esti.mode=2,cl.mode=0)
   majorcl(tr,shrink.mode=1,esti.mode=2)
   
3. Classification of simulation data by Dirichlet distribution (this simulation model tries to mimic mussel data, therefore this function will not work correctly after download due to error on missing "pomus" data. For research convinience, please refer to to the function code)
   da<-simotu.diri40(num,p=693,t=3,ref=c(1,214,490,512,513),seed=1234,unif.min=0.9,unif.max=0.95)  ##num is sample size for each target
   tr<-trotu(da,Target %in% c("target1","target2","target3"),thr=0,target="Target",
         pairs=c("target1","target2","target3"),del.otu=F,del.sam=T,nvar=75,
         lambda=seq(0.001,0.3,by=0.01),nsim=50,seed=F,nfold=5, nsampling=20, 
         test.per=0.2,norm.mode=1,shrink.mode=1,esti.mode=2,cl.mode=0)
   majorcl(tr,shrink.mode=1,esti.mode=2)
