
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

#prepare data
ind<- c("green", "wet")
# xmat<- data.matrix(path[,ind])
# npix<- path$n
# y<- path$dt


#model args
ngibbs=5000
# nburn=ngibbs/2
# w=0.01
# MaxIter=1000

#priors
var.betas=rep(10, length(ind))


#Run model
set.seed(1234)
plan(multisession)

progressr::with_progress({
  res<- resist(data = path, covs = ind, priors = var.betas, ngibbs = ngibbs)
})
future:::ClusterRegistry("stop")  #close all threads and memory used
# takes 55 sec to run (for 5000 iter)




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
betas<- data.frame(id = paste(rep(names(store.betas), each = length(var.betas)),
                              rep(ind, length(unique(path$id)))),
                   betas)
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
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  # scale_color_manual(values = pal) +
  scale_color_fish_d(option = "Scarus_tricolor") +
  theme_bw() +
  coord_flip() +
  labs(x="", y="") +
  theme(axis.text = element_text(size = 14),
        panel.grid = element_blank())

# ggsave("Giant Armadillo Resistance_coeffs.png", width = 9, height = 5,
#        units = "in", dpi = 300)



######################
### Export results ###
######################


# write.csv(post, "Giant Armadillo Resistance Results.csv", row.names = F)
