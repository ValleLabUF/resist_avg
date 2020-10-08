rm(list=ls())
set.seed(3)

#generate covariate from ar1 process
nobs=100000
nsteps=200 #how many times each covariate 
cov1=rep(runif(nobs/nsteps,min=-1,max=1),each=nsteps)
cov2=rep(runif(nobs/nsteps,min=-1,max=1),each=nsteps)
cov3=rep(runif(nobs/nsteps,min=-1,max=1),each=nsteps)
xmat=cbind(1,cov1,cov2,cov3)
nomes=paste0('x',1:3)
colnames(xmat)[-1]=nomes

#parameters 
betas=betas.true=c(0.5,0.1,-0.3,0)
b.gamma.true=b.gamma=1.3

#simulate data
media=exp(xmat%*%betas)
a1=media*b.gamma
delta.t=rgamma(nobs,a1,b.gamma)
hist(delta.t)

time.tot=cumsum(delta.t)
fim=as.data.frame(xmat[,nomes])
fim$time.tot=time.tot
fim$min1=NA
fim$npix=NA
fim$select=0

#select data every 30 seconds
target.time=30
tmp=abs(fim$time.tot-target.time)
ind1=0
fim$seg.id=NA
xavg=matrix(NA,10000,ncol(xmat)-1)
oo=1
for (i in 1:100000){
  print(i)
  
  #get new ind0 and ind1
  ind0=ind1+1
  max1=min(c(ind0+1000,nobs))
  tmp=abs(fim$time.tot[ind0:max1]-target.time*i)
  min1=min(tmp)
  ind1=which(tmp==min1)+(ind0-1)

  #use this information to set seg.id and select
  seq1=ind0:ind1
  fim$seg.id[seq1]=i
  fim$select[ind1]=1
  fim$min1[ind1]=min1
  fim$n[ind1]=length(seq1)

  #get covariates
  xavg[oo,]=apply(fim[seq1,nomes],2,mean); oo=oo+1
  
  #stop if finished
  if (ind1==nobs) break
}
hist(fim$min1)
table(fim$n)
sum(fim$select)

colnames(xavg)=nomes
plot(xmat[,'x1'],xavg[fim$seg.id,'x1'])
plot(xmat[,'x2'],xavg[fim$seg.id,'x2'])

fim$dt=fim$select*target.time
cond=fim$select==1
fim1=fim[cond,]

setwd('U:\\GIT_models\\resist_avg')
cond=!is.na(xavg[,1])
fim2=cbind(fim1[,c('n','dt')],xavg[cond,])
write.csv(fim2,'fake data resistance model.csv',row.names=F)