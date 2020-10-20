library(Rcpp)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(tidyr)

source('gibbs_resist_avg.R')
source('gibbs_resist_avg_func.R')
source('slice_b_gamma.R')
source('slice_betas.R')


# N and S IDs separated
path.N<- read.csv('N Armadillo Resistance Data_LULC.csv', as.is=T)
path.S<- read.csv('S Armadillo Resistance Data_LULC.csv', as.is=T)

path.N$dt<- path.N$dt/60  #convert to min from sec
path.S$dt<- path.S$dt/60


# Filter data for only steps with 3 >= dt >= 7 min
path.N<- path.N[path.N$dt >= 3 & path.N$dt <= 7 & !is.na(path.N$dt),]
path.S<- path.S[path.S$dt >= 3 & path.S$dt <= 7 & !is.na(path.S$dt),]


# Filter data by behavior (foraging or transit)
path.N.forage<- path.N %>% 
  filter(state == "Foraging")
path.N.transit<- path.N %>% 
  filter(state == "Transit")

path.S.forage<- path.S %>% 
  filter(state == "Foraging")
path.S.transit<- path.S %>% 
  filter(state == "Transit")



# Center and Scale covariates 
path.N.forage<- path.N.forage %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)
path.N.transit<- path.N.transit %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)

path.S.forage<- path.S.forage %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)
path.S.transit<- path.S.transit %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)






######################
### North Pantanal ###
######################

## Foraging

names(path.N.forage)[2:7]<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
xmat<- data.matrix(path.N.forage[,ind])
npix<- path.N.forage$n


#get y soma
y=path.N.forage$dt

  
#model args
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=rep(10,ncol(xmat))


#Run model
set.seed(2)
mod.forage_N<- gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                     w=w,MaxIter=MaxIter,npix=npix)
# takes 8 min to run (for 1000 iter)





## Transit

names(path.N.transit)[2:7]<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
xmat<- data.matrix(path.N.transit[,ind])
npix<- path.N.transit$n


#get y soma
y=path.N.transit$dt


#model args
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=rep(10,ncol(xmat))


#Run model
set.seed(2)
mod.transit_N<- gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                            w=w,MaxIter=MaxIter,npix=npix)
# takes 8 min to run (for 1000 iter)









######################
### South Pantanal ###
######################

# ### Group 'Campo' and 'Pasto' together to try fixing autocorr issue
# path.S.forage<- path.S.forage %>%
#   mutate_at("lulc", ~recode(., '5' = '1'))
# path.S.transit<- path.S.transit %>%
#   mutate_at("lulc", ~recode(., '5' = '1'))



names(path.S.forage)[2:6]<- c("Field", "Forest", "Water", "Pasture", "Road")
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S.forage[,ind])
npix<- path.S.forage$n


#get y soma
y=path.S.forage$dt


#model args
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=rep(10,ncol(xmat))


#Run model
set.seed(2)
mod.forage_S<- gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                            w=w,MaxIter=MaxIter,npix=npix)
# takes 5 min to run (for 1000 iter)







## Transit

names(path.S.transit)[2:6]<- c("Field", "Forest", "Water", "Pasture", "Road")
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S.transit[,ind])
npix<- path.S.transit$n


#get y soma
y=path.S.transit$dt


#model args
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=rep(10,ncol(xmat))


#Run model
set.seed(2)
mod.transit_S<- gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                            w=w,MaxIter=MaxIter,npix=npix)
# takes 1.5 min to run (for 1000 iter)





###################################################
### Check Posterior and Perform Model Selection ###
###################################################

######################
### North Pantanal ###
######################

## Foraging

#store results
store.llk.forage_N<- mod.forage_N$llk
store.b.forage_N<- mod.forage_N$b.gamma
store.betas.forage_N<- mod.forage_N$betas

#look at overall convergence
plot(store.llk.forage_N, type='l')
abline(v=nburn, col='red')
plot(store.llk.forage_N[(nburn + 1):ngibbs], type='l')
acf(store.llk.forage_N[(nburn + 1):ngibbs])

plot(store.b.forage_N, type='l')
plot(store.b.forage_N[(nburn + 1):ngibbs], type='l')
acf(store.b.forage_N[(nburn + 1):ngibbs])

#look at convergence betas
par(mfrow=c(2,2))
nbetas.forage_N<- ncol(mod.forage_N$betas)
for (i in 1:nbetas.forage_N){
  plot(mod.forage_N$betas[,i], type='l')  
}

for (i in 1:nbetas.forage_N){
  plot(mod.forage_N$betas[(nburn + 1):ngibbs, i], type='l')  
}
par(mfrow=c(1,1),mar=rep(3,4))

## Traceplots all indicate that convergence has been reached






# AIC_mcmc = function(llk, npar) {
#   (-2 * llk) + (2*npar) 
# }
# 
# AIC_N_NoInt<- AIC_mcmc(llk = mean(store.llk_N[(nburn+1):ngibbs,]), npar = 6)
# AIC_N_Int<- AIC_mcmc(llk = mean(store.llk_N2[(nburn+1):ngibbs,]), npar = 7)
# 
# AIC_N_NoInt - AIC_N_Int  #model 2 (w/ interaction) is much better






## Transit

#store results
store.llk.transit_N<- mod.transit_N$llk
store.b.transit_N<- mod.transit_N$b.gamma
store.betas.transit_N<- mod.transit_N$betas

#look at overall convergence
plot(store.llk.transit_N, type='l')
abline(v=nburn, col='red')
plot(store.llk.transit_N[(nburn + 1):ngibbs], type='l')
acf(store.llk.transit_N[(nburn + 1):ngibbs])

plot(store.b.transit_N, type='l')
plot(store.b.transit_N[(nburn + 1):ngibbs], type='l')
acf(store.b.transit_N[(nburn + 1):ngibbs])

#look at convergence betas
par(mfrow=c(2,2))
nbetas.transit_N<- ncol(mod.transit_N$betas)
for (i in 1:nbetas.transit_N){
  plot(mod.transit_N$betas[,i], type='l')  
}

for (i in 1:nbetas.transit_N){
  plot(mod.transit_N$betas[(nburn + 1):ngibbs, i], type='l')  
}
par(mfrow=c(1,1),mar=rep(3,4))

## Traceplots all indicate that convergence has been reached










######################
### South Pantanal ###
######################

## Foraging

#store results
store.llk.forage_S<- mod.forage_S$llk
store.b.forage_S<- mod.forage_S$b.gamma
store.betas.forage_S<- mod.forage_S$betas

#look at overall convergence
plot(store.llk.forage_S, type='l')
abline(v=nburn, col='red')
plot(store.llk.forage_S[(nburn + 1):ngibbs], type='l')
acf(store.llk.forage_S[(nburn + 1):ngibbs])

plot(store.b.forage_S, type='l')
plot(store.b.forage_S[(nburn + 1):ngibbs], type='l')
acf(store.b.forage_S[(nburn + 1):ngibbs])

#look at convergence betas
par(mfrow=c(2,2))
nbetas.forage_S<- ncol(mod.forage_S$betas)
for (i in 1:nbetas.forage_S){
  plot(mod.forage_S$betas[,i], type='l')  
}

for (i in 1:nbetas.forage_S){
  plot(mod.forage_S$betas[(nburn + 1):ngibbs, i], type='l')  
}
par(mfrow=c(1,1),mar=rep(3,4))

## Traceplots all indicate that convergence has been reached







## Transit

#store results
store.llk.transit_S<- mod.transit_S$llk
store.b.transit_S<- mod.transit_S$b.gamma
store.betas.transit_S<- mod.transit_S$betas

#look at overall convergence
plot(store.llk.transit_S, type='l')
abline(v=nburn, col='red')
plot(store.llk.transit_S[(nburn + 1):ngibbs], type='l')
acf(store.llk.transit_S[(nburn + 1):ngibbs])

plot(store.b.transit_S, type='l')
plot(store.b.transit_S[(nburn + 1):ngibbs], type='l')
acf(store.b.transit_S[(nburn + 1):ngibbs])

#look at convergence betas
par(mfrow=c(2,2))
nbetas.transit_S<- ncol(mod.transit_S$betas)
for (i in 1:nbetas.transit_S){
  plot(mod.transit_S$betas[,i], type='l')  
}

for (i in 1:nbetas.transit_S){
  plot(mod.transit_S$betas[(nburn + 1):ngibbs, i], type='l')  
}
par(mfrow=c(1,1),mar=rep(3,4))

## Traceplots all indicate that convergence has been reached








### FIND WAY TO EXPORT AND SAVE RESULTS TO USE IN ANALYSES AND DATA VIZ

# Export results

# write.csv(, "N Armadillo Resistance Results_dispersal.csv", row.names = F)