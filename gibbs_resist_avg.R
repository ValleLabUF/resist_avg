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
  
  
  #progress bar
  # pb<- progress::progress_bar$new(
  #   format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
  #   total = ngibbs, clear = FALSE, width = 100)

  
  for (i in 1:ngibbs){
    # pb$tick()  #create progress bar
    
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
#--------------------------------------------


### create wrapper for gibbs sampler

resist = function(data, covs, priors, ngibbs) {
  ## data  A data frame containing all pertinent vars; required vars include 'id', 'n' and 'dt',
  ##       although 'id' can go by any name
  ## covs  A vector of the column names that are storing the relevant covariates
  ## priors  A vector of the priors on the variance of the beta coeffs
  ## ngibbs integer. The number of iterations of the Gibbs sampler
  
  
  dat.list<- bayesmove::df_to_list(dat = data, ind = "id")
  
  #prepare data
  xmat<- purrr::map(dat.list, ~.x[,c("id",covs)])
  npix<- purrr::map(dat.list, ~purrr::pluck(.x, "n"))
  y<- purrr::map(dat.list, ~purrr::pluck(.x, "dt"))
  
  #tuning params
  w=0.01
  MaxIter=1000
  nburn = ngibbs/2
  
  #set up progress bar
  p<- progressr::progressor(steps = length(xmat))
  
  #run model
  tictoc::tic()
  dat.res<- furrr::future_pmap(list(y,xmat,npix),
    function(a, b, c) {
      tmp<- gibbs_resist(y=a, xmat=data.matrix(b[,-1]), ngibbs=ngibbs, nburn=nburn,
                         var.betas=priors, w=w, MaxIter=MaxIter, npix=c)
      p()  #for progress bar
      
      tmp
      },
    .options = furrr_options(seed = TRUE))
  tictoc::toc()
  
  
  dat.res  
}
