library(R2jags)
library(tictoc)
library(MCMCvis)


set.seed(2)

#get data
dat=read.csv('fake data resistance model.csv',as.is=T)


#define the model
model = function(){
  #likelihood
  for (i in 1:nobs){
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

#prepare data for jags
nobs=nrow(dat)
dt=dat$dt
n=dat$n
xmat=cbind(1,dat[,c('x1','x2','x3')])
dat1=list(nobs=nobs,dt=dt,n=n,xmat=xmat,nparam=ncol(xmat))

#set parameters to track
params=c('betas','b')
    
# MCMC settings 
n.iter <- 1000
n.thin <- 10  
n.burnin <- n.iter/2
n.chains <- 3


    
#fit the model
tic()
res = jags(model.file = model, parameters.to.save = params, data = dat1,
               n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter,
               n.thin = n.thin, DIC = TRUE)
toc()
# takes 2 min to run

res

MCMCsummary(res)
MCMCtrace(res, ind = TRUE, iter = 50, pdf = FALSE)
MCMCplot(res, excl = "deviance")
