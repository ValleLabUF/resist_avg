library(R2jags)
library(tictoc)
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyr)



# N and S IDs separated
path.N<- read.csv('N Armadillo Resistance Data_LULC.csv', as.is=T)
path.S<- read.csv('S Armadillo Resistance Data_LULC.csv', as.is=T)

path.N$dt<- path.N$dt/60  #convert to min from sec
path.S$dt<- path.S$dt/60


# Filter data for only steps with 3 >= dt >= 7 min
path.N<- path.N[path.N$dt >= 3 & path.N$dt <= 7 & !is.na(path.N$dt),]
path.S<- path.S[path.S$dt >= 3 & path.S$dt <= 7 & !is.na(path.S$dt),]



# Remove Burrow behavior observations
path.N<- path.N %>% 
  filter(state != 'Burrow')
path.S<- path.S %>% 
  filter(state != 'Burrow')



# Center and Scale covariates 
path.N<- path.N %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)

path.S<- path.S %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)






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
path.N.list<- bayesmove::df_to_list(path.N, "id")
path.N<- path.N.list[order(sapply(path.N.list, nrow), decreasing=TRUE)] %>% #reorder IDs
  bind_rows()
id<- as.numeric(factor(path.N$id, levels = unique(path.N$id)))
n.id<- max(id)
nobs<- nrow(path.N)
dt<- path.N$dt
n<- path.N$n
names(path.N)[2:7]<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
xmat<- data.matrix(path.N[,ind])
dat1<- list(nobs=nobs, dt=dt, n=n, xmat=xmat, nparam=ncol(xmat), id=id, n.id=n.id)


#set parameters to track
params=c('betas','b','b0')


#MCMC settings 
n.iter <- 2000
n.thin <- 5  
n.burnin <- n.iter/2
n.chains <- 3


#Run model
set.seed(2)

tic()
res.N = jags(model.file = model, parameters.to.save = params, data = dat1,
           n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter,
           n.thin = n.thin, DIC = TRUE)
toc()
# takes 36 min to run 2000 iterations


res.N

MCMCsummary(res.N)
MCMCtrace(res.N, ind = TRUE, iter = 100, pdf = FALSE)
par(mfrow=c(1,1))
MCMCplot(res.N, excl = "deviance")






######################
### South Pantanal ###
######################


#prepare data for jags
path.S.list<- bayesmove::df_to_list(path.S, "id")
path.S<- path.S.list[order(sapply(path.S.list, nrow), decreasing=TRUE)] %>% #reorder IDs
  bind_rows()
id<- as.numeric(factor(path.S$id, levels = unique(path.S$id)))
n.id<- max(id)
nobs<- nrow(path.S)
dt<- path.S$dt
n<- path.S$n
names(path.S)[2:6]<- c("Field", "Forest", "Water", "Pasture", "Road")
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S[,ind])
dat1<- list(nobs=nobs, dt=dt, n=n, xmat=xmat, nparam=ncol(xmat), id=id, n.id=n.id)


#set parameters to track
params=c('betas','b','b0')


#MCMC settings 
n.iter <- 2000
n.thin <- 5  
n.burnin <- n.iter/2
n.chains <- 3


#Run model
set.seed(2)

tic()
res.S = jags(model.file = model, parameters.to.save = params, data = dat1,
           n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter,
           n.thin = n.thin, DIC = TRUE)
toc()
# takes 19 min to run 2000 iterations


res.S

MCMCsummary(res.S)
MCMCtrace(res.S, ind = TRUE, iter = 100, pdf = FALSE)
par(mfrow=c(1,1))
MCMCplot(res.S, excl = "deviance")













### FIND WAY TO EXPORT AND SAVE RESULTS TO USE IN ANALYSES AND DATA VIZ

# Export results
res.N.summ<- res.N$BUGSoutput$summary
res.S.summ<- res.S$BUGSoutput$summary

#JAGS model of foraging/transit behaviors
# write.csv(res.N.summ, "N Armadillo Resistance Results.csv", row.names = F)
# write.csv(res.S.summ, "S Armadillo Resistance Results.csv", row.names = F)
