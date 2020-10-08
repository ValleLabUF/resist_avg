Doubling_Betas=function(xmat,betas1,w,b.gamma,yslice,MaxIter,y,target.p,
                        sd.betas,npix){
  betasLo=betasHi=betas1
  betasLo[target.p]=betasLo[target.p]-w*runif(1)
  betasHi[target.p]=betasLo[target.p]+w
  
  #calculate llk
  ylo=sum(get.llk(betas=betasLo,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix))+
      sum(dnorm(betasLo,mean=0,sd=sd.betas,log=T))
  yhi=sum(get.llk(betas=betasHi,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix))+
      sum(dnorm(betasHi,mean=0,sd=sd.betas,log=T))
  cond=F
  
  #keep doubling until ylo<yslice and yhi<yslice
  oo=0
  while ((ylo>yslice) & (oo<MaxIter)){
    betasLo[target.p]=betasLo[target.p]-w
    ylo=sum(get.llk(betas=betasLo,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix))+
        sum(dnorm(betasLo,mean=0,sd=sd.betas,log=T))
    oo=oo+1
  }
  if (oo >= MaxIter) cond=T
  oo=0
  while ((yhi>yslice) & (oo<MaxIter)){
    betasHi[target.p]=betasHi[target.p]+w
    yhi=sum(get.llk(betas=betasHi,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix))+
        sum(dnorm(betasHi,mean=0,sd=sd.betas,log=T))
    oo=oo+1
  }
  if (oo >= MaxIter) cond=T
  if (cond) return(c(betas1[target.p],betas1[target.p]))
  
  c(betasLo[target.p],betasHi[target.p])
}
#-------------------------------------------
Shrink_Sample_Betas=function(rango1,yslice,MaxIter,betas1,y,xmat,
                             b.gamma,target.p,sd.betas,npix){
  diff1=rango1[2]-rango1[1]
  betas1.orig=betas1
  yfim=-Inf
  oo=0
  while ((yfim<yslice) & (diff1 > 0.00000001) & (oo<MaxIter)){
    x=runif(1,min=rango1[1],max=rango1[2])
    betas1[target.p]=x
    yfim=sum(get.llk(betas=betas1,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix))+
         sum(dnorm(betas1,mean=0,sd=sd.betas,log=T))
    if (yfim<yslice){ #shrink the slice if x falls outside
      DistLo=abs(rango1[1]-x)
      DistHi=abs(rango1[2]-x)
      if (DistLo>DistHi) rango1[1]=x
      if (DistLo<DistHi) rango1[2]=x
      diff1=rango1[2]-rango1[1]
    }
    oo=oo+1
  }
  if (oo >= MaxIter | diff1 <= 0.00000001) return(betas1.orig)
  betas1
}
#-------------------------------------------
Sample_betas=function(nparam,xmat,y,betas,b.gamma,sd.betas,w,MaxIter,npix){
  for (j in 1:nparam){
     #define upper bound
    upper1=sum(get.llk(betas=betas,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix))+
           sum(dnorm(betas,mean=0,sd=sd.betas,log=T))
    yslice=upper1-rexp(1);
      
    #define slice
    rango1=Doubling_Betas(xmat=xmat,betas1=betas,w=w,b.gamma=b.gamma,
                          yslice=yslice,MaxIter=MaxIter,y=y,
                          target.p=j,sd.betas=sd.betas,npix=npix)
      
    #sample this particular parameter
    betas=Shrink_Sample_Betas(rango1=rango1,yslice=yslice,MaxIter=MaxIter,
                               betas1=betas,y=y,xmat=xmat,b.gamma=b.gamma,
                               target.p=j,sd.betas=sd.betas,npix=npix) 
  }
  betas
}