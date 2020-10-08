store.llk=mod.res$llk
store.b=mod.res$b.gamma
store.betas=mod.res$betas

#look at correlation
k=cor(cbind(store.b,store.betas))
k[k < 0.5 & k > -0.5]=NA
diag(k)=NA
k

#look at overall convergence
par(mfrow=c(1,1))
plot(store.llk,type='l')
nburn=700
abline(v=nburn,col='red')
plot(store.llk[nburn:ngibbs],type='l')

#look at convergence betas
par(mfrow=c(2,2),mar=rep(2,4))
for (i in 1:3) plot(store.betas[1:ngibbs,i],type='l')
for (i in 1:3) {
  plot(store.betas[nburn:ngibbs,i],type='l')
  abline(h=betas.true[i],col='red')
}

#look at convergence of b.gammas
par(mfrow=c(1,1))
plot(store.b,type='l')
abline(h=b.gamma.true,col='red')

#look at betas
par(mfrow=c(1,1),mar=rep(3,4))
betas.estim=store.betas[ngibbs,]
rango=range(c(betas.estim,betas.true))
plot(betas.true,betas.estim,xlim=rango,ylim=rango)
lines(rango,rango,col='red')