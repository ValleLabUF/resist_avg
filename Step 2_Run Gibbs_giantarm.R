
library(R2jags)
library(tictoc)
library(MCMCvis)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(furrr)
library(purrr)
library(coda)
library(stringr)
library(fishualize)


# Load data
path<- read.csv('Giant Armadillo Resistance Data.csv', as.is=T)
path$dt<- path$dt/60  #convert to min from sec
path$date<- as_datetime(path$date)


# Filter data for only steps with 5 >= dt >= 9 min
path<- path[path$dt >= 5 & path$dt <= 9 & !is.na(path$dt),]


# Remove state col and rows where is.na(greenness)
path<- path %>% 
  # dplyr::select(-state) %>% 
  filter(!is.na(green))


# Center and scale covariates 
# path$ndwi<- as.numeric(scale(path$ndwi, center = TRUE, scale = TRUE))
path<- path %>% 
  mutate(across(elev:wet, scale)) %>% 
  dplyr::select(-elev)  #remove elev to see if it removes weird predictions


# Add cols for season
# month1<- month.abb[month(path$date)]
# month1<- factor(month1, levels = month.abb[c(5:12,1)])
# season1<- ifelse(month1 %in% c(month.abb[3:5]), "Fall",
#                     ifelse(month1 %in% c(month.abb[6:8]), "Winter",
#                            ifelse(month1 %in% c(month.abb[9:11]), "Spring", "Summer")))
# season1<- factor(season1, levels = c("Fall","Winter","Spring","Summer"))
# table(season1)  #since Winter has largest N, this will be reference
# 
# season.mat<- model.matrix(~season1 + 0)
# season.mat<- season.mat[,-2]
# colnames(season.mat)<- c("Fall","Spring","Summer")
# 
# path2<- cbind(path, season.mat)



# Add cols for interaction between season and NDVI
# path3<- path2 %>% 
#   mutate(Fall.ndvi = ndvi*Fall,
#          Spring.ndvi = ndvi*Spring,
#          Summer.ndvi = ndvi*Summer)





#################
### Run model ###
#################

#define the model
model = function(){
  #likelihood
  for (i in 1:nobs){
    mu[i] <- n[i]*exp(b0[id[i]] + inprod(xmat[i,],betas[id[i], 1:nparam]))
    a[i] <- mu[i]*b
    dt[i] ~ dgamma(a[i],b)
  } 
  
  
  #priors
  b ~ dunif(0,1000)
  
  for (k in 1:nparam){  #random effect on slopes
    for (l in 1:n.id) {  
      betas[l,k] ~ dnorm(mu_slope, tau2_slope)
    }
  }
  mu_slope ~ dnorm(0,0.1)
  tau2_slope ~ dgamma(0.1,0.1)
  
  for (j in 2:n.id){  #random effect on intercept
    b0[j] ~ dnorm(0, tau2_int)
  }
  b0[1]<- 0  #ref ID is set to 0
  # mu.id ~ dnorm(0,0.1)
  tau2_int ~ dgamma(0.1,0.1)
}





#prepare data for jags
path.list<- bayesmove::df_to_list(path, "id")
path2<- path.list[order(sapply(path.list, nrow), decreasing=TRUE)] %>% #reorder IDs
  bind_rows()
id<- as.numeric(factor(path2$id, levels = unique(path2$id)))
n.id<- max(id)
nobs<- nrow(path2)
dt<- path2$dt
n<- path2$n
ind<- c("green","wet")
xmat<- data.matrix(path2[,ind])
dat1<- list(nobs=nobs, dt=dt, n=n, xmat=xmat, nparam=ncol(xmat), id=id, n.id=n.id)


#set parameters to track
params=c('betas','b','b0')


#MCMC settings 
n.iter <- 5000
n.thin <- 5  
n.burnin <- n.iter/2
n.chains <- 3


#Run model
set.seed(123)

tic()
res = jags(model.file = model, parameters.to.save = params, data = dat1,
             n.chains = n.chains, n.burnin = n.burnin, n.iter = n.iter,
             n.thin = n.thin, DIC = TRUE)
toc()
# takes 14.5 min to run 5000 iterations


res

MCMCsummary(res)
MCMCtrace(res, ind = TRUE, iter = 500, pdf = FALSE)
par(mfrow=c(1,1))
MCMCplot(res, excl = "deviance")

res.summ<- res$BUGSoutput$summary




### Make (pretty) caterpillar plot
betas<- data.frame(res.summ[2:22,c(1,3,7)])
betas$coeff<- c(rep('int', n.id),
                rep('green', n.id),
                rep('wet', n.id))
betas$id<- rep(unique(path2$id), 3)
names(betas)[2:3]<- c("lower", "upper")
betas$coeff<- factor(betas$coeff, levels = c("int","green","wet"))


ggplot(data=betas, aes(x=coeff, y=mean, ymin=lower, ymax=upper, color = id)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  scale_x_discrete(labels = c("Intercept","Greenness", "Wetness")) +
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

id1<- unique(path2$id)


## Greenness
green.res<- list()

for (i in 1:length(unique(path2$id))) {
  
  #Generate sequence along green
  rango1<- path2 %>% 
    filter(id == id1[i]) %>% 
    dplyr::select(green) %>% 
    range()
  seq.green<- seq(rango1[1], rango1[2], length.out = 100)
  
  
  #Create design matrix where 0s added for all other vars besides green
  design.mat<- cbind(betas[betas$id == id1[i] & betas$coeff == "int", "mean"],
                     seq.green,
                     0)
  
  # Take cross-product of design matrix with betas and exponentiate to calc response
  ind<- which(id1[i] == betas$id)  #indicator by ID
  
  y.mu<- exp(design.mat %*% betas$mean[ind])
  y.low<- exp(design.mat %*% betas$lower[ind])
  y.up<- exp(design.mat %*% betas$upper[ind])
  
  
  # Add results to data frame
  y.mu.df<- data.frame(x = seq.green,
                       y = y.mu,
                       ymin = y.low,
                       ymax = y.up,
                       id = id1[i])
  
  green.res[[i]]<- y.mu.df
}

green.res.df<- bind_rows(green.res)

# Plot relationship
ggplot(data = green.res.df) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, fill = id), alpha =  0.3) +
  geom_line(aes(x, y, color = id), size = 1) +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  labs(x = "\nStandardized Greenness", y = "Time Spent per Cell (min)\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))





## Wetness
wet.res<- list()

for (i in 1:length(unique(path2$id))) {
  
  #Generate sequence along wet
  rango1<- path2 %>% 
    filter(id == id1[i]) %>% 
    dplyr::select(wet) %>% 
    range()
  seq.wet<- seq(rango1[1], rango1[2], length.out = 100)
  
  
  #Create design matrix where 0s added for all other vars besides wet
  design.mat<- cbind(betas[betas$id == id1[i] & betas$coeff == "int", "mean"],
                     0,
                     seq.wet)
  
  # Take cross-product of design matrix with betas and exponentiate to calc response
  ind<- which(id1[i] == betas$id)  #indicator by ID
  
  y.mu<- exp(design.mat %*% betas$mean[ind])
  y.low<- exp(design.mat %*% betas$lower[ind])
  y.up<- exp(design.mat %*% betas$upper[ind])
  
  
  # Add results to data frame
  y.mu.df<- data.frame(x = seq.wet,
                       y = y.mu,
                       ymin = y.low,
                       ymax = y.up,
                       id = id1[i])
  
  wet.res[[i]]<- y.mu.df
}

wet.res.df<- bind_rows(wet.res)

# Plot relationship
ggplot(data = wet.res.df) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, fill = id), alpha =  0.3) +
  geom_line(aes(x, y, color = id), size = 1) +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  labs(x = "\nStandardized Wetness", y = "Time Spent per Cell (min)\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))




######################
### Export results ###
######################


# write.csv(betas, "Giant Armadillo Resistance Results.csv", row.names = F)
