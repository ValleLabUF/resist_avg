
library(R2jags)
library(tictoc)
library(MCMCvis)
library(tidyverse)
library(lubridate)
library(furrr)
library(coda)
library(fishualize)


# Load data
path<- read.csv('Giant Armadillo Resistance Data.csv', as.is=T)
path$dt<- path$dt/60  #convert to min from sec
path$date<- as_datetime(path$date)


# Filter data for only steps with 5 >= dt >= 9 min
path<- path[path$dt >= 5 & path$dt <= 9 & !is.na(path$dt),]


# Remove state col and rows where is.na(evi) and dist == 0
path<- path %>% 
  # dplyr::select(-state) %>% 
  filter(!is.na(evi)) %>% 
  filter(dist > 0)


# Center and scale covariates 
# path$ndwi<- as.numeric(scale(path$ndwi, center = TRUE, scale = TRUE))
path<- path %>% 
  mutate(across(evi, scale))




#################
### Run model ###
#################

#define the model
model = function(){
  #likelihood
  for (i in 1:nobs){
    tmp[i] <- exp(betas[id[i],1] + betas[id[i],2] * evi[i])
    a[i] <- b[id[i]] * dist[i] * tmp[i] 
    dt[i] ~ dgamma(a[i], b[id[i]])
  }
  
  #random effects on betas (intercept and slope) and rate param
  rho ~ dunif(-1,1)
  sig[1] ~ dunif(0,5)
  sig[2] ~ dunif(0,5) 
  
  Sigma[1,1] <- pow(sig[1], 2) 
  Sigma[2,2] <- pow(sig[2], 2) 
  Sigma[1,2] <- sig[1] * sig[2] * rho
  Sigma[2,1] <- Sigma[1,2]
  mu <- c(0, 0)
  
  for (j in 1:n.id) {
    betas[j, 1:2] ~ dmnorm.vcov(mu[1:2], Sigma[1:2,1:2])
    b[j] ~ dunif(0,100)
  }
  
}





#prepare data for jags
# path.list<- bayesmove::df_to_list(path, "id")
# path2<- path.list[order(sapply(path.list, nrow), decreasing=TRUE)] %>% #reorder IDs
#   bind_rows()
id<- as.numeric(factor(path$id, levels = unique(path$id)))
n.id<- max(id)
nobs<- nrow(path)
dt<- path$dt
dist<- path$dist
evi<- as.numeric(path$evi)
# ind<- 'evi'
# xmat<- data.matrix(path2[,ind])
dat1<- list(nobs=nobs, dt=dt, dist=dist, evi=evi, id=id, n.id=n.id)


#set parameters to track
params=c('b','betas','sig','rho')


#MCMC settings 
# n.iter <- 5000
# n.thin <- 5  
# n.burnin <- n.iter/2
# n.chains <- 3


#Run model
# set.seed(1)

tic()
res = jags.parallel(model.file = model, parameters.to.save = params, data = dat1,
             n.chains = 3, n.burnin = 3000, n.iter = 6000,
             n.thin = 5, DIC = TRUE, jags.seed = 123)
toc()
# takes 4 min to run 5000 iterations


res

MCMCsummary(res)
MCMCtrace(res, ind = TRUE, iter = 600, pdf = FALSE)
par(mfrow=c(1,1))
MCMCplot(res, excl = "deviance")

res.summ<- res$BUGSoutput$summary




### Make (pretty) caterpillar plot
params<- data.frame(res.summ[c(1:21),c(1,3,7)])
params$coeff<- rep(c('b','int','evi'), each = n.id)
params$id<- rep(unique(path$id), 3)
params$id<- factor(params$id, levels = unique(params$id))
names(params)[2:3]<- c("lower", "upper")
params$coeff<- factor(params$coeff, levels = c('b',"int","evi"))


ggplot(data=params, aes(x=coeff, y=mean, ymin=lower, ymax=upper, color = id)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  # geom_errorbar(data = betas[22:23,], aes(x=coeff, y=mean, ymin=lower, ymax=upper),
  #               position = position_dodge(0.55), width = 0, size = 0.75, color = "black") +
  # geom_point(data = betas[22:23,], aes(x=coeff, y=mean, ymin=lower, ymax=upper),
  #            position = position_dodge(0.55), size=2, color = "black") +
  scale_x_discrete(labels = c("b","Intercept","EVI")) +
  scale_color_fish_d("", option = "Scarus_tricolor") +
  theme_bw() +
  coord_flip() +
  labs(x="", y="") +
  theme(axis.text = element_text(size = 14),
        panel.grid = element_blank())

# ggsave("Giant Armadillo Resistance_coeffs.png", width = 9, height = 5,
#        units = "in", dpi = 300)






###########################################
### Viz partial responses to covs by ID ###
###########################################

id1<- unique(path$id)
dist1<- mean(path$dist)
b1<- params[params$coeff == "b",]
betas<- params[params$coeff != "b",]

## EVI
evi.res<- list()

for (i in 1:length(unique(path$id))) {
  
  tmp<- path %>% 
    filter(id == id1[i])
  
  #Generate sequence along evi
  rango1<- tmp %>% 
    dplyr::select(evi) %>% 
    range()
  seq.evi<- seq(rango1[1], rango1[2], length.out = nrow(tmp))
  
  
  #Create design matrix where 0s added for all other vars besides green
  design.mat<- cbind(1, seq.evi)
  
  # Take cross-product of design matrix with betas and exponentiate to calc response
  ind<- which(id1[i] == betas$id)  #indicator by ID
  
  y.mu<- b1[i,1] * dist1 * exp(design.mat %*% betas$mean[ind])
  y.low<- b1[i,2] * dist1 * exp(design.mat %*% betas$lower[ind])
  y.up<- b1[i,3] * dist1 * exp(design.mat %*% betas$upper[ind])
  
  
  # Add results to data frame
  y.mu.df<- data.frame(x = seq.evi,
                       y = y.mu,
                       ymin = y.low,
                       ymax = y.up,
                       id = id1[i])
  
  evi.res[[i]]<- y.mu.df
}

evi.res.df<- bind_rows(evi.res)

# Plot relationship
ggplot(data = evi.res.df) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, fill = id), alpha =  0.3) +
  geom_line(aes(x, y, color = id), size = 1) +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  labs(x = "\nStandardized EVI", y = "Time Spent per Cell (min)\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))







######################
### Export results ###
######################


# write.csv(params, "Giant Armadillo Resistance Results.csv", row.names = F)
