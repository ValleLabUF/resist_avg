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


# # Filter data by behavior (foraging or transit)
# path.N.forage<- path.N %>% 
#   filter(state == "Foraging")
# path.N.transit<- path.N %>% 
#   filter(state == "Transit")
# 
# path.S.forage<- path.S %>% 
#   filter(state == "Foraging")
# path.S.transit<- path.S %>% 
#   filter(state == "Transit")


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

names(path.N)[2:7]<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
xmat<- data.matrix(path.N[,ind])
npix<- path.N$n


#get y soma
y=path.N$dt

  
#model args
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=rep(10,ncol(xmat))


#Run model
set.seed(2)
mod_N<- gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                     w=w,MaxIter=MaxIter,npix=npix)
# takes 11 min to run (for 1000 iter)









######################
### South Pantanal ###
######################


names(path.S)[2:6]<- c("Field", "Forest", "Water", "Pasture", "Road")
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S[,ind])
npix<- path.S$n


#get y soma
y=path.S$dt


#model args
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=rep(10,ncol(xmat))


#Run model
set.seed(2)
mod_S<- gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                            w=w,MaxIter=MaxIter,npix=npix)
# takes 5 min to run (for 1000 iter)









###################################################
### Check Posterior and Perform Model Selection ###
###################################################

######################
### North Pantanal ###
######################

#store results
store.llk_N<- mod_N$llk
store.b_N<- mod_N$b.gamma
store.betas_N<- mod_N$betas

#look at overall convergence
plot(store.llk_N, type='l')
abline(v=nburn, col='red')
plot(store.llk_N[(nburn + 1):ngibbs], type='l')
acf(store.llk_N[(nburn + 1):ngibbs])

plot(store.b_N, type='l')
plot(store.b_N[(nburn + 1):ngibbs], type='l')
acf(store.b_N[(nburn + 1):ngibbs])

#look at convergence betas
par(mfrow=c(2,2))
nbetas_N<- ncol(mod_N$betas)
for (i in 1:nbetas_N){
  plot(mod_N$betas[,i], type='l')  
}

for (i in 1:nbetas_N){
  plot(mod_N$betas[(nburn + 1):ngibbs, i], type='l')  
}
par(mfrow=c(1,1),mar=rep(3,4))

## Traceplots all indicate that convergence has been reached









######################
### South Pantanal ###
######################

#store results
store.llk_S<- mod_S$llk
store.b_S<- mod_S$b.gamma
store.betas_S<- mod_S$betas

#look at overall convergence
plot(store.llk_S, type='l')
abline(v=nburn, col='red')
plot(store.llk_S[(nburn + 1):ngibbs], type='l')
acf(store.llk_S[(nburn + 1):ngibbs])

plot(store.b_S, type='l')
plot(store.b_S[(nburn + 1):ngibbs], type='l')
acf(store.b_S[(nburn + 1):ngibbs])

#look at convergence betas
par(mfrow=c(2,2))
nbetas_S<- ncol(mod_S$betas)
for (i in 1:nbetas_S){
  plot(mod_S$betas[,i], type='l')  
}

for (i in 1:nbetas_S){
  plot(mod_S$betas[(nburn + 1):ngibbs, i], type='l')  
}
par(mfrow=c(1,1),mar=rep(3,4))

## Traceplots all indicate that convergence has been reached









### FIND WAY TO EXPORT AND SAVE RESULTS TO USE IN ANALYSES AND DATA VIZ

# Export results

# write.csv(, "N Armadillo Resistance Results_dispersal.csv", row.names = F)