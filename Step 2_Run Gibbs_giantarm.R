
library(R2jags)
library(tictoc)
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)



# Load data
path<- read.csv('Giant Armadillo Resistance Data.csv', as.is=T)
path$dt<- path$dt/60  #convert to min from sec
path$date<- as_datetime(path$date)


# Filter data for only steps with 5 >= dt >= 9 min
path<- path[path$dt >= 5 & path$dt <= 9 & !is.na(path$dt),]


# Remove state col and rows where is.na(NDVI)
path<- path %>% 
  dplyr::select(-state) %>% 
  filter(!is.na(ndvi))


# Center and scale covariates 
path$ndvi<- as.numeric(scale(path$ndvi, center = TRUE, scale = TRUE))




######################
### North Pantanal ###
######################


#define the model
model = function(){
  #likelihood
  for (i in 1:nobs){
    mu[i] <- n[i]*exp(b0[id[i]] + inprod(xmat[i,],betas[1:nparam]))
    a[i] <- mu[i]*b
    dt[i] ~ dgamma(a[i],b)
  } 
  
  
  #priors
  b ~ dunif(0,1000)
  
  for (k in 1:nparam){
    betas[k] ~ dnorm(0,0.1)
  }
  
  for (j in 2:n.id){
    b0[j] ~ dnorm(0, tau2)
  }
  b0[1]<- 0  #ref ID is set to 0
  # mu.id ~ dnorm(0,0.1)
  tau2 ~ dgamma(0.1,0.1)
}


#prepare data for jags
path.list<- bayesmove::df_to_list(path, "id")
path<- path.list[order(sapply(path.list, nrow), decreasing=TRUE)] %>% #reorder IDs
  bind_rows()
id<- as.numeric(factor(path$id, levels = unique(path$id)))
n.id<- max(id)
nobs<- nrow(path)
dt<- path$dt
n<- path$n
# names(path)[3:6]<- c("...", "...", "...", "...")
ind<- c("ndvi","X1", "X2", "X3", "X4")
xmat<- data.matrix(path[,ind])
dat1<- list(nobs=nobs, dt=dt, n=n, xmat=xmat, nparam=ncol(xmat), id=id, n.id=n.id)


#set parameters to track
params=c('betas','b','b0')


#MCMC settings 
n.iter <- 1000
n.thin <- 5  
n.burnin <- n.iter/2
n.chains <- 3


#Run model
set.seed(1234)

tic()
res = jags(model.file = model, parameters.to.save = params, data = dat1,
             n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter,
             n.thin = n.thin, DIC = TRUE)
toc()
# takes 15.5 min to run 1000 iterations


res

MCMCsummary(res)
MCMCtrace(res, ind = TRUE, iter = 100, pdf = FALSE)
par(mfrow=c(1,1))
MCMCplot(res, excl = "deviance")





######################
### Export results ###
######################

res.summ<- res$BUGSoutput$summary

# write.csv(res.summ, "Giant Armadillo Resistance Results.csv", row.names = F)
