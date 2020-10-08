gibbs_resist=function(y,xmat,ngibbs,nburn,var.betas,w,MaxIter,npix){
  n=nrow(xmat)  
  nparam=ncol(xmat)

  #initial parameters
  betas=matrix(0,nparam,1)
  betas[1]=mean(log(y/npix))
  b.gamma=1

  #priors
  sd.betas=sqrt(var.betas)

  #stuff for gibbs sampler
  store.betas=matrix(NA,ngibbs,nparam)
  store.b=matrix(NA,ngibbs,1)
  store.llk=matrix(NA,ngibbs,1)

  for (i in 1:ngibbs){
    print(i)
    # print(sum(abs(betas)))
    
    #sample betas
    betas=Sample_betas(nparam=nparam,xmat=xmat,y=y,betas=betas,
                       b.gamma=b.gamma,sd.betas=sd.betas,w=w,MaxIter=MaxIter,
                       npix=npix)
    # betas=betas.true
    
    #sample b.gamma
    b.gamma=Sample_bgamma(nparam=nparam,xmat=xmat,
                          y=y,betas=betas,b.gamma=b.gamma,
                          w=w,MaxIter=MaxIter,npix=npix)
    # b.gamma=1
    
    #get llk
    p=get.llk(betas=betas,xmat=xmat,y=y,b.gamma=b.gamma,npix=npix)
    llk1=sum(p)

    #store results
    store.betas[i,]=betas
    store.b[i]=b.gamma
    store.llk[i]=sum(llk1)
  }
  list(betas=store.betas,b.gamma=store.b,llk=store.llk)
}



