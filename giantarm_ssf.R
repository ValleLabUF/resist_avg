
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

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")
dat<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)

# Filter out observations where coords are NA
dat<- dat %>% 
  filter(!is.na(x))


dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))

dat<- dat %>% 
  rename(t = date) %>% 
  mutate_at("t", as_datetime)



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
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

#LULC
lulc<- raster('cheiann_UTM1.tif')
names(lulc)<- "lulc"

#NDWI
files<- list.files(getwd(), pattern = "*.grd$")
ndwi.filenames<- files[grep("ndwi", files)]
ndwi.season<- brick(ndwi.filenames[2])
ndwi.season<- resample(ndwi.season, lulc, method = "bilinear")
compareRaster(lulc, ndwi.season)

#set times (Z) for dynamic ndwi RasterBrick (still need to provide a value for lulc too)
ndwi<- setZ(ndwi.season, ymd(c("2019-03-01","2019-06-01","2019-09-01","2019-12-01")))

## will need to extract lulc and ndwi separately



##########################################
### Analysis of Giant Armadillo Tracks ###
##########################################

# Generate available steps for each ID and extract environ covars
dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% steps_by_burst() %>% 
                      random_steps(n_control = 30)))

dat_ssf2<- dat_ssf %>% 
  mutate(stps = map(stps, ~.x %>% 
                      extract_covariates_var_time(covariates = ndwi,
                                                              when = "after",
                                                              where = "end",
                                                              max_time = days(91),
                                                              name_covar = "ndwi") %>% 
                      extract_covariates(lulc))) %>% 
  select(id, stps) %>% unnest() %>% 
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste0(id, step_id_, sep = "-"))
dat_ssf2


# Remove obs where NDWI is NA
dat_ssf3<- dat_ssf2[-which(is.na(dat_ssf2$ndwi)),]


# Convert 'lulc' to a dummy variable
group<- factor(dat_ssf3$lulc)
lulc_dum<- model.matrix(~ group + 0)
colnames(lulc_dum)<- c("Forest", "Closed_Savanna", "Open_Savanna", "Floodable")

# Add lulc contrast to data frame
dat_ssf4<- cbind(dat_ssf3, lulc_dum)

# Center and scale NDWI
dat_ssf4$ndwi<- scale(dat_ssf4$ndwi)



### glmmTMB for fitting SSF

## First try fitting a model w/o random slopes (leave out 'Pasture' as reference level)
# TMBStruc.fix = glmmTMB(y ~ HQ + Fence + Water + Cane + Forest + ndvi_N +
#                          (1|step_id), 
#                        family=poisson, data=dat_ssf.N2, doFit=FALSE) 
# 
# # Then fix the standard deviation of the first random term, which is the `(1|step_id)` component  in the above model equation:
# TMBStruc.fix$parameters$theta[1] = log(1e3) 
# 
# # We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
# TMBStruc.fix$mapArg = list(theta=factor(c(NA)))
# 
# # Then fit the model and look at the results:
# glmm.TMB.fixed = glmmTMB:::fitTMB(TMBStruc.fix) 
# summary(glmm.TMB.fixed)



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




### Fit SSF using clogit()
ssf.mod<- dat_ssf4 %>% 
  fit_clogit(y ~ Closed_Savanna + Open_Savanna + Floodable + ndwi + strata(step_id) + 
               survival::cluster(id))
(soma<- summary(ssf.mod))



### Make caterpillar plot
soma.betas<- as.data.frame(soma$conf.int[-5,-2])
names(soma.betas)<- c("mean", "lower", "upper")
soma.betas$coeff<- rownames(soma.betas)
soma.betas$coeff<- factor(soma.betas$coeff, levels = soma.betas$coeff)


ggplot(data=soma.betas, aes(x=coeff, y=mean, ymin=lower, ymax=upper)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(position = position_dodge(0.75), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.75), size=2) +
  # scale_color_manual(values = pal) +
  # scale_color_fish_d(option = "Scarus_tricolor") +
  theme_bw() +
  coord_flip() +
  labs(x="", y="") +
  theme(axis.text = element_text(size = 14),
        panel.grid = element_blank())

# ggsave("Giant Armadillo Habitat Selection_coeffs.png", width = 6, height = 5,
#        units = "in", dpi = 300)



### Make spatial predictions from SSF

#extract beta coeffs (mean)
betas<- coef(ssf.mod)

#Center and scale NDWI raster
ndwi_s<- scale(ndwi, center = T, scale = T)


##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("ndwi","Forest", "Closed_Savanna", "Open_Savanna", "Floodable")
xmat<- dat_ssf4 %>% 
  filter(y == 1) %>% 
  dplyr::select(all_of(ind)) %>% 
  data.matrix()

lulc.mat<- model.matrix.lm(~factor(getValues(lulc)) + 0, na.action = "na.pass")
colnames(lulc.mat)<- ind[2:5]



# Combine LULC contrast w NDWI by season

ssfSurf<- vector("list", 4)
names(ssfSurf)<- names(ndwi_s)
for (i in 1:nlayers(ndwi_s)) {
  mat<- cbind(lulc.mat[,-1], ndwi = values(ndwi_s[[i]]))
  
  tmp<- lulc
  raster::values(tmp)<- exp(mat %*% betas[1:4])
  
  ssfSurf[[i]]<- tmp
}


# Create function for negative exponential transformation of habitat suitability to resistance via Keeley et al. (2016):
neg.exp.trans = function(rast, c) {
  tmp<- 100-99*((1-exp(-c*rast)) / (1-exp(-c)))
  tmp
}

# Scale SSF prob and take inverse for resistance
for (i in 1:length(ssfSurf)) {
  ssfSurf[[i]]<- ssfSurf[[i]]/ssfSurf[[i]]@data@max
}

# ssfSurfN_r<- neg.exp.trans(ssfSurfN_s, 8)


#create as data frame
ssfSurf.df<- map(ssfSurf, as.data.frame, xy=T) %>% 
  map(., ~rename(., sel = lulc)) %>% 
  bind_rows(., .id = "season")
ssfSurf.df$season<- factor(ssfSurf.df$season, levels = names(ssfSurf))

# ssfSurfN_r.df<- as.data.frame(ssfSurfN_r, xy=T)
# names(ssfSurfN_r.df)[3]<- "resist"



## Selection (by season)
ggplot() +
  geom_raster(data = ssfSurf.df, aes(x, y, fill = sel)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "black") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Habitat Selection") +
  theme_bw() +
  coord_equal() +
  facet_wrap(~ season, nrow = 1) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))

# ggsave("Giant Armadillo Habitat Selection_season.png", width = 9, height = 5,
#        units = "in", dpi = 300)



## calculate mean across seasons
ssfSurf.mean.df<- ssfSurf.df %>% 
  pivot_wider(., names_from = season, values_from = sel) %>% 
  mutate(., mu = rowMeans(select(., 3:6), na.rm = TRUE))

## Selection (average)
ggplot() +
  geom_raster(data = ssfSurf.mean.df, aes(x, y, fill = mu)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "black") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Habitat Selection") +
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

# ggsave("Giant Armadillo Habitat Selection.png", width = 9, height = 5,
#        units = "in", dpi = 300)



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


