
### Fit giant armadillo movement using SSF

library(glmmTMB)
library(lubridate)
# library(INLA)
library(tidyverse)
library(raster)
library(survival)
# library(TwoStepCLogit)
library(amt)
library(fishualize)


#################
### Load data ###
#################

set.seed(2021)

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")

dat<- read.csv('Giant Armadillo state estimates.csv', as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing) %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'), factor,
                levels = c("Slow-Turn","Slow-Unif","Exploratory","Transit","Unclassified"))
  )

dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))

dat<- dat %>% 
  rename(t = date) %>% 
  mutate_at("t", as_datetime, tz = "UTC")



# Modify data into format readable by 'amt'
dat_all <- dat %>% nest(-id) 
dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, x, y, t, crs = sp::CRS("+init=epsg:32721"))
  }))


# Summarize time steps
dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  select(id, sr) %>% 
  unnest

# Resample the tracks to 7 minutes with a tolerance of 2 minutes
dat1 <- dat_all %>% 
  mutate(dat_clean = map(trk, ~ {
    .x %>% track_resample(rate = minutes(7), tolerance = seconds(120))
  }))



#Load in environmental covariates

#Tasseled Cap Greenness
green<- brick('GiantArm_tcgreen_season.grd')
green_s<- green  #for scaled greenness
values(green_s)<- scale(values(green_s))

#Tasseled Cap Wetness
wet<- brick('GiantArm_tcwet_season.grd')
wet_s<- wet  #for scaled wetness
values(wet_s)<- scale(values(wet_s))
compareRaster(green_s, wet_s)


#set times (Z) for dynamic greenness and wetness RasterBricks
green_s<- setZ(green_s, ymd(c("2019-03-01","2019-06-01","2019-09-01","2019-12-01")))
wet_s<- setZ(wet_s, ymd(c("2019-03-01","2019-06-01","2019-09-01","2019-12-01")))

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")




##########################################
### Analysis of Giant Armadillo Tracks ###
##########################################

# Generate available steps for each ID and extract environ covars
dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps(n_control = 30)))

dat_ssf2<- dat_ssf %>% 
  mutate(stps = map(stps, ~.x %>% 
                      extract_covariates_var_time(covariates = green_s,
                                                              when = "after",
                                                              where = "end",
                                                              max_time = days(91),
                                                              name_covar = "green") %>% 
                      extract_covariates_var_time(covariates = wet_s,
                                                  when = "after",
                                                  where = "end",
                                                  max_time = days(91),
                                                  name_covar = "wet")
                    )) %>% 
  select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"))
dat_ssf2


# Remove obs where green is NA
dat_ssf3<- dat_ssf2[-which(is.na(dat_ssf2$green)),]



### glmmTMB for fitting SSF

# ssf.list<- vector("list", length(unique(dat_ssf4$id)))
# names(ssf.list)<- unique(dat$id)
# 
# for (i in 1:length(ssf.list)) {
#   print(i)
#   ## First try fitting a model w/o random slopes (leave out 'Pasture' as reference level)
#   TMBStruc.fix = glmmTMB(y ~ Closed_Savanna + Open_Savanna + Floodable + ndwi +
#                            (1|step_id),
#                          family=poisson,
#                          data=dat_ssf4[dat_ssf4$id == i,],
#                          doFit=FALSE)
#   
#   # Then fix the standard deviation of the first random term, which is the `(1|step_id)` component  in the above model equation:
#   TMBStruc.fix$parameters$theta[1] = log(1e2)
#   
#   # We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
#   TMBStruc.fix$mapArg = list(theta=factor(c(NA)))
#   
#   # Then fit the model and look at the results:
#   glmm.TMB.fixed = glmmTMB:::fitTMB(TMBStruc.fix)
#   
#   #wrangle and store results
#   soma <- summary(glmm.TMB.fixed)
#   tmp<- data.frame(soma$coefficients$cond[,c("Estimate","Std. Error")])
#   tmp$coef.names<- rownames(tmp)
#   rownames(tmp)<- 1:nrow(tmp)
#   ssf.list[[i]]<- tmp[-1,]
# }


# Now the same model using `glmmTMB()`. Note that we do not need an overall intercept in this model, because the stratum-specific intercepts are (almost) freely estimated due to the large, fixed variance. Again start to set up the model without fitting the model:
# TMBStruc <- glmmTMB(y ~ -1 + HQ + Fence + Water + Cane + Forest + ndvi_N + 
#                       (0 + ndvi_N | id) + 
#                       (0 + HQ | id) + (0 + Fence | id) + (0 + Water | id) +
#                       (0 + Cane | id) + (0 + Forest | id) +
#                       (1|step_id),
#                     family=poisson, data = dat_ssf.N2, doFit=FALSE) 
# 
# # Set the value of the standard deviation of the first random effect (here (1|step_id)):
# TMBStruc$parameters$theta[7] <- log(1e3) 
# 
# # Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
# TMBStruc$mapArg <- list(theta=factor(c(1:6,NA)))
# 
# # Fit the model and look at the summary:
# glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
# summary(glmm.TMB.random)

# 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
# confint(glmm.TMB.random)


# Note: It is currently a problem of glmmTMB that the confint() function only shows a table of the length equal to the number of parameters estimated. As the variance for str_ID was not estimated but fixed to 10^6, but is still listed, the last variance component (here the one for (0 + forest | id)) is not shown. This can be solved by moving the component (1|step_id) to the last position in the formula, and then replace the respective lines by 
# TMBStruc$parameters$theta[1] = log(1e3) 
# TMBStruc$mapArg = list(theta=factor(c(1, NA)))




### Fit SSF using clogit() by ID

ssf.list<- vector("list", length(unique(dat_ssf3$id)))
names(ssf.list)<- unique(dat$id)

for (i in 1:length(ssf.list)) {
  #fit model
  ssf.mod <- dat_ssf3[dat_ssf3$id == i,] %>%
    fit_clogit(y ~ green + wet + strata(step_id))
  
  #wrangle and store results
  soma <- summary(ssf.mod)
  tmp<- data.frame(soma$coefficients[,c("coef","se(coef)")])
  tmp$coef.names<- rownames(tmp)
  rownames(tmp)<- 1:nrow(tmp)
  ssf.list[[i]]<- tmp
}



### Make caterpillar plot
soma.betas<- bind_rows(ssf.list, .id = "id")
names(soma.betas)[2:3]<- c("coef","se.coef")
soma.betas<- soma.betas %>% 
  mutate(lower = coef - (1.96 * se.coef),
         upper = coef + (1.96 * se.coef),
         .before = coef.names)
# names(soma.betas)<- c("mean", "lower", "upper")
# soma.betas$coeff<- rownames(soma.betas)
soma.betas$coef.names<- factor(soma.betas$coef.names, levels = unique(soma.betas$coef.names))


ggplot(data=soma.betas, aes(x=coef.names, y=coef, ymin=lower, ymax=upper, color = id)) +
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

# ggsave("Giant Armadillo Habitat Selection_coeffs.png", width = 9, height = 5,
#        units = "in", dpi = 300)



### Make spatial predictions from SSF

#extract beta coeffs (mean)
betas<- soma.betas$coef


##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("green","wet")
xmat<- dat_ssf3 %>% 
  filter(y == 1) %>% 
  dplyr::select(all_of(ind)) %>% 
  data.matrix()



# ## currently ssf.mod is fit for  a single ID; last stored model is for 'tex'
# 
# lulc.mat2<- data.frame(lulc.mat)[,-1]
# covars<- cbind(lulc.mat2, ndwi = values(ndwi_s$Fall), step_id = NA)
# 
# pred.fall<- predict(object = ssf.mod$model, newdata = covars, type = "risk", reference = "sample")
# pred.fall.prob<- pred.fall / (1 + pred.fall)



### Calculate Resistance Surface by ID ###

ssfSurf<- vector("list", 4)
names(ssfSurf)<- names(green_s)
id1<- unique(soma.betas$id)

for (j in 1:nlayers(green_s)) {
  print(names(green_s)[j])
  
  tmp<- list()
  
  cov.mat<- cbind(green = raster::values(green_s[[j]]), wet = raster::values(wet_s[[j]]))
  
  for (i in 1:length(id1)) {
    print(i)
    
    ssf.res<- green[[1]]
    w.hat<- exp(cov.mat %*% soma.betas[which(soma.betas$id == id1[i]), "coef"])
    raster::values(ssf.res)<- w.hat/(1 + w.hat)
    # resistSurf<- resistSurf * 60  #convert from min to sec
    
    #create as data frame
    ssf.res.df<- as.data.frame(ssf.res, xy=T) %>%
      mutate(id = i)
    names(ssf.res.df)[3]<- "sel"
    
    tmp[[i]]<- ssf.res.df
  }
  
  names(tmp)<- id1
  ssfSurf[[j]]<- tmp
  
}


ssfSurf.df<- map(ssfSurf, ~bind_rows(., .id = "id")) %>%
  bind_rows(., .id = "season")
ssfSurf.df$season<- factor(ssfSurf.df$season, levels = unique(ssfSurf.df$season))








# # Combine LULC contrast w NDWI by season
# 
# ssfSurf<- vector("list", 4)
# names(ssfSurf)<- names(ndwi_s)
# for (i in 1:nlayers(ndwi_s)) {
#   mat<- cbind(lulc.mat[,-1], ndwi = values(ndwi_s[[i]]))
#   
#   tmp<- lulc
#   raster::values(tmp)<- exp(mat %*% betas[1:4])
#   
#   ssfSurf[[i]]<- tmp
# }


# Create function for negative exponential transformation of habitat suitability to resistance via Keeley et al. (2016):
neg.exp.trans = function(rast, c) {
  tmp<- 100-99*((1-exp(-c*rast)) / (1-exp(-c)))
  tmp
}

# # Scale SSF prob and take inverse for resistance
# for (i in 1:length(ssfSurf)) {
#   ssfSurf[[i]]<- ssfSurf[[i]]/ssfSurf[[i]]@data@max
# }

# ssfSurfN_r<- neg.exp.trans(ssfSurfN_s, 8)


#create as data frame
# ssfSurf.df<- map(ssfSurf, as.data.frame, xy=T) %>% 
#   map(., ~rename(., sel = lulc)) %>% 
#   bind_rows(., .id = "season")
# ssfSurf.df$season<- factor(ssfSurf.df$season, levels = names(ssfSurf))

# ssfSurfN_r.df<- as.data.frame(ssfSurfN_r, xy=T)
# names(ssfSurfN_r.df)[3]<- "resist"



## Selection (by season and ID)
ggplot() +
  geom_raster(data = ssfSurf.df, aes(x, y, fill = sel)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "chartreuse") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Habitat Selection") +
  theme_bw() +
  coord_equal() +
  facet_grid(season ~ id) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))

# ggsave("Giant Armadillo Habitat Selection_IDxSeason_facet.png", width = 9, height = 7.5,
#        units = "in", dpi = 300)


#remove predictions for seasons in which ID has no observations; will only calculate mean habitat selection based on seasons w/ data
# seasons<- names(ndwi)
# id.tmp<- list()
# for (i in 1:length(unique(ssfSurf.df$id))) {
#   cond<- unique(dat[dat$id == unique(dat$id)[i], "season"])
#   ind1<- which(seasons %in% cond)
#   
#   id.tmp[[i]]<- ssfSurf.df %>% 
#     filter(id == unique(ssfSurf.df$id)[i]) %>% 
#     filter(season == seasons[ind1])
# }
# ssfSurf.df2<- bind_rows(id.tmp)


tmp<- bayesmove::df_to_list(dat, "id") %>% 
  purrr::map(., . %>% 
               dplyr::select(season) %>% 
               # unique() %>% 
               unlist())

weights<- tmp %>% 
  map(., ~{table(.x)/length(.x)}) %>% 
  map(., ~{.x[.x !=0]})


#calculate mean resistance across seasons
ssf.mean.id<- ssfSurf.df %>% 
  bayesmove::df_to_list(., "id") %>% 
  map2(., tmp, ~{.x %>% 
         filter(season %in% unique(.y))}) %>% 
  map(., bayesmove::df_to_list, "season") %>% 
  map_depth(., 2, pluck, "sel")

ssf.mean.id2<- ssf.mean.id %>% 
  map(., bind_cols) %>% 
  map2(., weights, ~apply(.x, 1, function(x) weighted.mean(x, .y, na.rm = TRUE))) %>% 
  unlist()

ssf.mean.id3<- cbind(ssfSurf$Fall$blanca[,c("x","y")], sel = ssf.mean.id2) %>% 
  mutate(id = rep(names(ssf.mean.id), each = ncell(green)), .before = "x")


## Mean across seasons (where observations recorded)
ggplot() +
  geom_raster(data = ssf.mean.id3, aes(x, y, fill = sel)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "chartreuse") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Habitat Selection") +
  theme_bw() +
  coord_equal() +
  facet_wrap( ~ id, nrow = 1) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))

# ggsave("Giant Armadillo Habitat Selection_ID_facet.png", width = 8.5, height = 4,
#        units = "in", dpi = 300)



# ## calculate mean across seasons
# ssfSurf.mean.df<- ssfSurf.df %>% 
#   pivot_wider(., names_from = season, values_from = sel) %>% 
#   mutate(., mu = rowMeans(select(., 3:6), na.rm = TRUE))

#calculate mean and variance across IDs
ssf.pop<- ssf.mean.id3 %>% 
  bayesmove::df_to_list(., "id") %>% 
  map_depth(., 1, pluck, "sel")

weights2<- table(dat$id)/nrow(dat)  #don't need to reorder for this script
weights2_red<- table(dat$id)[c(1:2,5:7)]/nrow(dat %>% filter(id != 'gala' & id != 'mafalda'))

ssf.pop.mean<- ssf.pop %>% 
  bind_cols() %>% 
  apply(., 1, function(x)  weighted.mean(x, weights2, na.rm = TRUE))

ssf.pop.mean_red<- ssf.pop %>%   #w/ removal of gala and mafalda
  bind_cols() %>% 
  dplyr::select(-c(gala, mafalda)) %>% 
  apply(., 1, function(x)  weighted.mean(x, weights2_red, na.rm = TRUE))

ssf.pop.var<- ssf.pop %>% 
  bind_cols() %>% 
  apply(., 1, var, na.rm = TRUE)

ssf.pop.var_red<- ssf.pop %>% 
  bind_cols() %>% 
  dplyr::select(-c(gala, mafalda)) %>% 
  apply(., 1, var, na.rm = TRUE)


ssf.pop2<- cbind(ssfSurf$Fall$blanca[,c("x","y")], mu = ssf.pop.mean,
                    sig = ssf.pop.var)

ssf.pop2_red<- cbind(ssfSurf$Fall$blanca[,c("x","y")], mu = ssf.pop.mean_red,
                 sig = ssf.pop.var_red)


## Selection (average; all IDs)
ggplot() +
  geom_raster(data = ssf.pop2, aes(x, y, fill = mu)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "black") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Mean of Population") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))

# ggsave("Giant Armadillo Habitat Selection_mean.png", width = 9, height = 5,
#        units = "in", dpi = 300)

## Selection (average; reduced IDs)
ggplot() +
  geom_raster(data = ssf.pop2_red, aes(x, y, fill = mu)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "black") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Mean of Population") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))


## Selection (variance; all IDs)
ggplot() +
  geom_raster(data = ssf.pop2, aes(x, y, fill = sig)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "chartreuse") +
  scale_fill_viridis_c("Selection", option = "viridis",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Variance of Population") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))

# ggsave("Giant Armadillo Habitat Selection_var.png", width = 9, height = 5,
#        units = "in", dpi = 300)


## Selection (variance; reduced IDs)
ggplot() +
  geom_raster(data = ssf.pop2_red, aes(x, y, fill = sig)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "chartreuse") +
  scale_fill_viridis_c("Selection", option = "viridis",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Variance of Population") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))



## Resistance (negative exponential transformation; c=8)
# ggplot() +
#   geom_raster(data = ssfSurfN_r.df, aes(x, y, fill = resist)) +
#   geom_path(data = dat %>% filter(region == "N"),
#             aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
#   scale_fill_viridis_c("Resistance", option = "inferno",
#                        na.value = "transparent") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   labs(x="Easting", y="Northing", title = "North Pantanal Resistance") +
#   theme_bw() +
#   coord_equal() +
#   theme(legend.position = "bottom",
#         axis.title = element_text(size = 18),
#         axis.text = element_text(size = 10),
#         strip.text = element_text(size = 16, face = "bold"),
#         plot.title = element_text(size = 22),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12)) +
#   guides(fill = guide_colourbar(barwidth = 30, barheight = 1))

# ggsave("N Pantanal Resistance.png", width = 9, height = 5, units = "in", dpi = 300)




### Export summary results

# write.csv(ssf.pop2_red, "Giant Armadillo SSF summary results.csv", row.names = F)
# write.csv(soma.betas, "SSF coeffs.csv", row.names = F)
