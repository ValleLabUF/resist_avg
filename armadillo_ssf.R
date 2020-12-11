
### Fit 3-banded armadillo movement using SSF

library(glmmTMB)
library(lubridate)
# library(INLA)
library(tidyverse)
library(raster)
# library(survival)
# library(TwoStepCLogit)
library(amt)


#################
### Load data ###
#################

dat<- read.csv("Armadillo HMM Results.csv", header = T, sep = ",")
dat<- dat %>% 
  rename(id = ID, t = date) %>% 
  mutate_at("id", as.character) %>% 
  mutate_at("state", as.factor) %>% 
  mutate_at("state", ~recode(., '1' = "Burrow",
                             '2' = "Foraging", '3' = "Transit"))
dat$t<- as_datetime(dat$t)

# Filter for only the active movements (foraging/transit)
dat2<- dat %>% filter(state == "Foraging" | state == "Transit")


# Modify data into format readable by 'amt'
dat_all <- dat2 %>% nest(-c(id, region)) 
dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, x, y, t, crs = sp::CRS("+init=epsg:32721"))
  }))


# Summarize time steps
dat_all %>% mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  select(id, sr) %>% 
  unnest

# Resample the tracks to 5 minutes with a tolerance of 1 minute
dat1 <- dat_all %>% 
  mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = minutes(5), tolerance = seconds(60))
}))

# Separate by region
dat.N<- dat1 %>% filter(region == "N")
dat.S<- dat1 %>% filter(region == "S")


#Load in environmental covariates
lulcN<- raster('lulc_N.tif')
ndviN<- raster('ndvi_N.tif')
covars.N<- stack(lulcN, ndviN)

lulcS<- raster('lulc_S.tif')
ndviS<- raster('ndvi_S.tif')
covars.S<- stack(lulcS, ndviS)




#####################################
### Analysis of N Pantanal Tracks ###
#####################################

# Generate available steps for each ID and extract environ covars
dat_ssf.N <- dat.N %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps() %>% extract_covariates(covars.N))) %>% 
  select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"))
dat_ssf.N


# Remove obs where NDVI (and LULC) are NA
dat_ssf.N<- dat_ssf.N[-which(is.na(dat_ssf.N$ndvi_N)),]


# Convert 'lulc_N' to a dummy variable
group<- factor(dat_ssf.N$lulc_N)
lulc_N_dum<- model.matrix(~ group + 0)
colnames(lulc_N_dum)<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")

# Add lulc_N contrast to data frame
dat_ssf.N2<- cbind(dat_ssf.N, lulc_N_dum)

# Center and scale NDVI and step length
dat_ssf.N2$ndvi_N<- scale(dat_ssf.N2$ndvi_N)



### glmmTMB for fitting SSF

## First try fitting a model w/o random slopes (leave out 'Pasture' as reference level)
TMBStruc.fix = glmmTMB(y ~ HQ + Fence + Water + Cane + Forest + ndvi_N +
                         (1|step_id), 
                       family=poisson, data=dat_ssf.N2, doFit=FALSE) 

# Then fix the standard deviation of the first random term, which is the `(1|step_id)` component  in the above model equation:
TMBStruc.fix$parameters$theta[1] = log(1e3) 

# We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
TMBStruc.fix$mapArg = list(theta=factor(c(NA)))

# Then fit the model and look at the results:
glmm.TMB.fixed = glmmTMB:::fitTMB(TMBStruc.fix) 
summary(glmm.TMB.fixed)



# Now the same model using `glmmTMB()`. Note that we do not need an overall intercept in this model, because the stratum-specific intercepts are (almost) freely estimated due to the large, fixed variance. Again start to set up the model without fitting the model:
TMBStruc <- glmmTMB(y ~ -1 + HQ + Fence + Water + Cane + Forest + ndvi_N + 
                      (0 + ndvi_N | id) + 
                      (0 + HQ | id) + (0 + Fence | id) + (0 + Water | id) +
                      (0 + Cane | id) + (0 + Forest | id) +
                      (1|step_id),
                    family=poisson, data = dat_ssf.N2, doFit=FALSE) 

# Set the value of the standard deviation of the first random effect (here (1|step_id)):
TMBStruc$parameters$theta[7] <- log(1e3) 

# Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBStruc$mapArg <- list(theta=factor(c(1:6,NA)))

# Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)

# 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)


# Note: It is currently a problem of glmmTMB that the confint() function only shows a table of the length equal to the number of parameters estimated. As the variance for str_ID was not estimated but fixed to 10^6, but is still listed, the last variance component (here the one for (0 + forest | id)) is not shown. This can be solved by moving the component (1|step_id) to the last position in the formula, and then replace the respective lines by 
TMBStruc$parameters$theta[1] = log(1e3) 
TMBStruc$mapArg = list(theta=factor(c(1, NA)))




### Fit SSF using clogit() only on "transit" state
ssf.modN<- dat_ssf.N2 %>% 
  fit_clogit(y ~ HQ + Fence + Water + Cane + Forest + ndvi_N + strata(step_id) + 
               survival::cluster(id))
summary(ssf.modN)




### Make spatial predictions from SSF

#extract beta coeffs (mean)
betas_N<- coef(ssf.modN)

#Center and scale NDVI raster
ndviN_s<- scale(ndviN, center = T, scale = T)


##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "ndvi_N")
xmat<- dat_ssf.N2 %>% 
  filter(y == 1) %>% 
  dplyr::select(all_of(ind)) %>% 
  data.matrix()

lulcN.mat<- model.matrix(~factor(getValues(lulcN)) + 0)
colnames(lulcN.mat)<- ind[1:6]

  

# Combine LULC contrast w NDVI
N.mat<- cbind(lulcN.mat[,-1], ndvi = values(ndviN_s))

ssfSurfN<- lulcN
raster::values(ssfSurfN)<- exp(N.mat %*% betas_N[1:6])


# Create function for negative exponential transformation of habitat suitability to resistance via Keeley et al. (2016):
neg.exp.trans = function(rast, c) {
  tmp<- 100-99*((1-exp(-c*rast)) / (1-exp(-c)))
  tmp
}

# Scale SSF prob and take inverse for resistance
ssfSurfN_s<- ssfSurfN/ssfSurfN@data@max
ssfSurfN_r<- neg.exp.trans(ssfSurfN_s, 8)


#create as data frame
ssfSurfN_s.df<- as.data.frame(ssfSurfN_s, xy=T)
names(ssfSurfN_s.df)[3]<- "sel"

ssfSurfN_r.df<- as.data.frame(ssfSurfN_r, xy=T)
names(ssfSurfN_r.df)[3]<- "resist"



## Selection
ggplot() +
  geom_raster(data = ssfSurfN_s.df, aes(x, y, fill = sel)) +
  geom_path(data = dat %>% filter(region == "N"),
            aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "North Pantanal Habitat Selection") +
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

setwd("~/Documents/Snail Kite Project/Local Presentation/")
# ggsave("N Pantanal Habitat Selection.png", width = 9, height = 5, units = "in", dpi = 300)


## Resistance (negative exponential transformation; c=8)
ggplot() +
  geom_raster(data = ssfSurfN_r.df, aes(x, y, fill = resist)) +
  geom_path(data = dat %>% filter(region == "N"),
            aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Resistance", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "North Pantanal Resistance") +
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

# ggsave("N Pantanal Resistance.png", width = 9, height = 5, units = "in", dpi = 300)







#####################################
### Analysis of S Pantanal Tracks ###
#####################################

# Generate available steps for each ID and extract environ covars
dat_ssf.S <- dat.S %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps() %>% extract_covariates(covars.S))) %>% 
  select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"))
dat_ssf.S


# Filter for only "transit" behavior

# Remove obs where "available" points have missing NDVI

dat_ssf.S2<- dat_ssf.S[-which(is.na(dat_ssf.S$ndvi_S)),]

# Convert 'lulc_S' to a dummy variable
group<- factor(dat_ssf.S2$lulc_S)
lulc_S_dum<- model.matrix(~ group + 0)
colnames(lulc_S_dum)<- c("Field", "Forest", "Water", "Pasture", "Road")

# Add lulc_S contrast to data frame
dat_ssf.S2<- cbind(dat_ssf.S2, lulc_S_dum)

# Center and scale NDVI
dat_ssf.S2$ndvi_S<- scale(dat_ssf.S2$ndvi_S)



### glmmTMB for fitting SSF

## First try fitting a model w/o random slopes (leave out 'Pasture' as reference level)
TMBStruc.fix = glmmTMB(y ~ Forest + Water + Pasture + Road + ndvi_S +
                         (1|step_id), 
                       family=poisson, data=dat_ssf.S2, doFit=FALSE) 

# Then fix the standard deviation of the first random term, which is the `(1|step_id)` component  in the above model equation:
TMBStruc.fix$parameters$theta[1] = log(1e3) 

# We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
TMBStruc.fix$mapArg = list(theta=factor(c(NA)))

# Then fit the model and look at the results:
glmm.TMB.fixed = glmmTMB:::fitTMB(TMBStruc.fix) 
summary(glmm.TMB.fixed)



# Now the same model using `glmmTMB()`. Note that we do not need an overall intercept in this model, because the stratum-specific intercepts are (almost) freely estimated due to the large, fixed variance. Again start to set up the model without fitting the model:
TMBStruc <- glmmTMB(y ~ -1 + HQ + Fence + Water + Cane + Forest + ndvi_S + 
                      (0 + ndvi_S | id) + 
                      (0 + HQ | id) + (0 + Fence | id) + (0 + Water | id) +
                      (0 + Cane | id) + (0 + Forest | id) +
                      (1|step_id),
                    family=poisson, data = dat_ssf.S2, doFit=FALSE) 

# Set the value of the standard deviation of the first random effect (here (1|step_id)):
TMBStruc$parameters$theta[7] <- log(1e3) 

# Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
TMBStruc$mapArg <- list(theta=factor(c(1:6,NA)))

# Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)

# 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)


# Note: It is currently a problem of glmmTMB that the confint() function only shows a table of the length equal to the number of parameters estimated. As the variance for str_ID was not estimated but fixed to 10^6, but is still listed, the last variance component (here the one for (0 + forest | id)) is not shown. This can be solved by moving the component (1|step_id) to the last position in the formula, and then replace the respective lines by 
TMBStruc$parameters$theta[1] = log(1e3) 
TMBStruc$mapArg = list(theta=factor(c(1, NA)))




### Fit SSF using clogit()
ssf.modS<- dat_ssf.S2 %>% 
  fit_clogit(y ~ Forest + Water + Pasture + Road + ndvi_S + strata(step_id) + 
               survival::cluster(id))
summary(ssf.modS)




### Make spatial predictions from SSF

#extract beta coeffs (mean)
betas_S<- coef(ssf.modS)

#Center and scale NDVI raster
ndviS_s<- scale(ndviS, center = T, scale = T)


##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "ndvi_S")
xmat<- dat_ssf.S2 %>% 
  filter(y == 1) %>% 
  dplyr::select(all_of(ind)) %>% 
  data.matrix()

lulcS.mat<- model.matrix(~factor(getValues(lulcS)) + 0)
colnames(lulcS.mat)<- ind[1:5]



# Combine LULC contrast w NDVI
S.mat<- cbind(lulcS.mat[,-1], ndvi = values(ndviS_s))

ssfSurfS<- lulcS
raster::values(ssfSurfS)<- exp(S.mat %*% betas_S[1:5])

# Scale SSF prob and take inverse for resistance
ssfSurfS_s<- ssfSurfS/ssfSurfS@data@max
ssfSurfS_r<- neg.exp.trans(ssfSurfS_s, c=8)


#create as data frame
ssfSurfS_s.df<- as.data.frame(ssfSurfS_s, xy=T)
names(ssfSurfS_s.df)[3]<- "sel"

ssfSurfS_r.df<- as.data.frame(ssfSurfS_r, xy=T)
names(ssfSurfS_r.df)[3]<- "resist"



## Selection
ggplot() +
  geom_raster(data = ssfSurfS_s.df, aes(x, y, fill = sel)) +
  geom_path(data = dat %>% filter(region == "S"),
            aes(x, y, group = id), alpha = 0.5, color = "blue") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Habitat Selection") +
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

# ggsave("S Pantanal Habitat Selection.png", width = 9, height = 5, units = "in", dpi = 300)


## Resistance (negative exponential transformation; c=8)
ggplot() +
  geom_raster(data = ssfSurfS_r.df, aes(x, y, fill = resist)) +
  geom_path(data = dat %>% filter(region == "S"),
            aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Resistance", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Resistance") +
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

# ggsave("S Pantanal Resistance.png", width = 9, height = 5, units = "in", dpi = 300)
