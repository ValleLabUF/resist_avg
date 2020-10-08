rm(list=ls())
library('Rcpp')
library('mvtnorm')
set.seed(92)

setwd('U:\\GIT_models\\resist_avg')
source('gibbs_resist_avg.R')
source('gibbs_resist_avg_func.R')
source('slice_b_gamma.R')
source('slice_betas.R')
dat=read.csv('fake data resistance model.csv',as.is=T)
ind=grep('x',colnames(dat))
xmat=data.matrix(cbind(1,dat[,ind])) #notice inclusion of intercept
npix=dat$n

#get y soma
y=dat$dt
ngibbs=1000
nburn=ngibbs/2
w=0.01
MaxIter=1000

#priors
var.betas=c(100,rep(10,ncol(xmat)-1))

mod.res=gibbs_resist(y=y,xmat=xmat,ngibbs=ngibbs,nburn=nburn,var.betas=var.betas,
                     w=w,MaxIter=MaxIter,npix=npix)