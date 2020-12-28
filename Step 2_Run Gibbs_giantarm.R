
# library(R2jags)
library(tictoc)
# library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(furrr)
library(purrr)
library(Rcpp)
library(mvtnorm)
library(coda)
library(stringr)
library(fishualize)

source('gibbs_resist_avg.R')
source('gibbs_resist_avg_func.R')
source('slice_b_gamma.R')
source('slice_betas.R')


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




#################
### Run model ###
#################

#prepare data
# names(path)[3:6]<- c("...", "...", "...", "...")
ind<- c("ndvi","X1", "X2", "X3", "X4")
# xmat<- data.matrix(path[,ind])
# npix<- path$n
# y<- path$dt


#model args
ngibbs=1000
# nburn=ngibbs/2
# w=0.01
# MaxIter=1000

#priors
var.betas=rep(10, length(ind))


#Run model
set.seed(1234)
plan(multisession)
res<- resist(data = path, covs = ind, priors = var.betas, ngibbs = ngibbs)
# takes 1.5 min to run (for 1000 iter)




#store results
store.llk<- map(res, ~pluck(., "llk"))
store.b<- map(res, ~pluck(., "b.gamma"))
store.betas<- map(res, ~pluck(., "betas"))

#look at overall convergence
llk<- t(bind_rows(store.llk))
llk<- data.frame(id = names(store.llk), llk)
bayesmove::traceplot(data = llk, ngibbs = ngibbs, type = "LML")
for (i in 1:length(store.llk)) {
  acf(store.llk[[i]][((ngibbs/2) + 1):ngibbs])
}

#look at convergence of b hyperparameter
b<- t(bind_rows(store.b))
b<- data.frame(id = names(store.b), b)
bayesmove::traceplot(data = b, ngibbs = ngibbs, type = "LML")
for (i in 1:length(store.b)) {
  acf(store.b[[i]][((ngibbs/2) + 1):ngibbs])
}

#look at convergence of betas
betas<- t(bind_rows(store.betas))
betas<- data.frame(id = paste(rep(names(store.betas), each = 5), rep(ind, 5)), betas)
bayesmove::traceplot(data = betas, ngibbs = ngibbs, type = "LML")
betas2<- map(store.betas, data.frame) %>% 
  bind_cols()
for (i in 1:ncol(betas2)) {
  acf(betas2[((ngibbs/2) + 1):ngibbs, i])
}



### Make caterpillar plot
betas3<- as.mcmc(betas2[((ngibbs/2)+1):ngibbs,])
hpd<- HPDinterval(betas3)
x.bar<- colMeans(betas3)
post<- data.frame(var = betas$id, mean = x.bar, hpd)
post<- cbind(str_split_fixed(string = post$var, pattern = " ", n = 2), post) %>% 
  dplyr::select(-var)
names(post)[1:2]<- c("id", "coeff")
post$coeff<- factor(post$coeff, levels = ind)


# pal<- as.character(palette.colors(length(unique(post$id))))

ggplot(data=post, aes(x=coeff, y=mean, ymin=lower, ymax=upper, color = id)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(position = position_dodge(0.75), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.75), size=2) +
  # scale_color_manual(values = pal) +
  scale_color_fish_d(option = "Scarus_tricolor") +
  theme_bw() +
  coord_flip() +
  labs(x="", y="") +
  theme(axis.text = element_text(size = 14),
        panel.grid = element_blank())



######################
### Export results ###
######################


# write.csv(post, "Giant Armadillo Resistance Results.csv", row.names = F)
