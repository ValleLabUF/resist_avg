Doubling_bgamma=function(media,w,b.gamma,yslice,MaxIter,y){
  b.gammaLo=b.gamma-w*runif(1)
  b.gammaLo=ifelse(b.gammaLo<0.0000001,0.0000001,b.gammaLo)
  b.gammaHi=b.gammaLo+w
  
  #calculate llk
  ylo=llk_bgamma(media=media,y=y,b.gamma=b.gammaLo)
  yhi=llk_bgamma(media=media,y=y,b.gamma=b.gammaHi)
  cond=F
  
  #keep doubling until ylo<yslice and yhi<yslice
  oo=0
  while ((ylo>yslice) & (oo<MaxIter)){
    b.gammaLo=b.gammaLo-w
    if (b.gammaLo<0.0000001) { #avoid negative values
      b.gammaLo=0.0000001
      break;
    }
    ylo=llk_bgamma(media=media,y=y,b.gamma=b.gammaLo)
    oo=oo+1
  }
  if (oo >= MaxIter) cond=T
  
  oo=0
  while ((yhi>yslice) & (oo<MaxIter)){
    b.gammaHi=b.gammaHi+w
    yhi=llk_bgamma(media=media,y=y,b.gamma=b.gammaHi)
    oo=oo+1
  }
  if (oo >= MaxIter) cond=T
  if (cond) return(c(b.gamma,b.gamma))
  
  c(b.gammaLo,b.gammaHi)
}
#-------------------------------------------
llk_bgamma=function(media,y,b.gamma){
  a1=b.gamma*media
  sum(dgamma(y,a1,b.gamma,log=T))
}
#-------------------------------------------
Shrink_Sample_bgamma=function(rango1,yslice,MaxIter,media,y,b.gamma){
  diff1=rango1[2]-rango1[1]
  yfim=-Inf
  oo=0
  while ((yfim<yslice) & (diff1 > 0.00000001) & (oo<MaxIter)){
    x=runif(1,min=rango1[1],max=rango1[2])
    yfim=llk_bgamma(media=media,y=y,b.gamma=x)
    if (yfim<yslice){ #shrink the slice if x falls outside
      DistLo=abs(rango1[1]-x)
      DistHi=abs(rango1[2]-x)
      if (DistLo>DistHi) rango1[1]=x
      if (DistLo<DistHi) rango1[2]=x
      diff1=rango1[2]-rango1[1]
    }
    oo=oo+1
  }
  if (diff1 <= 0.00000001 | oo >=MaxIter) return(b.gamma)
  x
}
#-------------------------------------------
Sample_bgamma=function(nparam,xmat,y,betas,b.gamma,w,MaxIter,npix){
  media1=npix*exp(xmat%*%betas)

  #define upper bound
  upper1=llk_bgamma(media=media1,y=y,b.gamma=b.gamma)
  yslice=upper1-rexp(1);
    
  #define slice
  rango1=Doubling_bgamma(media=media1,w=w,b.gamma=b.gamma,
                         yslice=yslice,MaxIter=MaxIter,y=y)

  #sample this particular parameter
  b.gamma=Shrink_Sample_bgamma(rango1=rango1,yslice=yslice,MaxIter=MaxIter,
                               media=media1,y=y,b.gamma=b.gamma)
  b.gamma
}