get.llk=function(betas,xmat,y,b.gamma,npix){
  media=npix*exp(xmat%*%betas)
  a.gamma1=b.gamma*media
  dgamma(y,a.gamma1,b.gamma,log=T)
}