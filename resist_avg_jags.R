model{
  for(i in 1:nobs){
    mu[i] <- n[i]*exp(inprod(xmat[i,],betas[1:nparam]))
    a[i] <- mu[i]*b
    dt[i] ~ dgamma(a[i],b)
  } 

  #priors
  b ~ dunif(0,1000)
  for (k in 1:nparam){
    betas[k] ~ dnorm(0,0.1)
  }
}
