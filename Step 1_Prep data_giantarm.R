### Load armadillo telemetry data

library(tidyverse)
library(sf)
library(raster)
library(lubridate)
library(sp)
library(furrr)
library(future)

source('helper functions.R')

#############################
### Import armadillo data ###
#############################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")
dat<- read.csv('Giant Armadillo state estimates.csv', as.is = T)
dat$date<- as_datetime(dat$date, tz = "UTC")
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



### Only importing environ covars that were used in the analysis w/ behavioral states (elevation, greenness, wetness)


#######################################
### Import Environmental Covariates ###
#######################################

#Tasseled Cap Greenness

green<- brick('GiantArm_tcgreen_season.grd')
green.df<- as.data.frame(green, xy = T)
green.df2<- pivot_longer(green.df, cols = -c(x,y), names_to = "season", values_to = "green")
green.df2$season<- factor(green.df2$season, levels = names(green))

ggplot() +
  geom_raster(data = green.df2, aes(x, y, fill = green)) +
  scale_fill_distiller("Greenness", palette = "Greens", na.value = "transparent",
                       limits = c(-2000,5000), direction = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold")) +
  facet_wrap(~ season)


#Tasseled Cap Wetness

wet<- brick('GiantArm_tcwet_season.grd')
compareRaster(green, wet)

wet.df<- as.data.frame(wet, xy = T)
wet.df2<- pivot_longer(wet.df, cols = -c(x,y), names_to = "season", values_to = "wet")
wet.df2$season<- factor(wet.df2$season, levels = names(wet))

#map TC Wetness
ggplot() +
  geom_raster(data = wet.df2, aes(x, y, fill = wet)) +
  scale_fill_distiller("Wetness", palette = "Blues", na.value = "transparent",
                       limits = c(-5000,250), direction = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold")) +
  facet_wrap(~ season)



#Elevation
dem<- raster('giantarm_dem.tif')
names(dem)<- 'elev'
dem<- resample(dem, green, method = "bilinear")
compareRaster(green, dem)
dem.df<- as.data.frame(dem, xy = TRUE)

ggplot() +
  geom_raster(data = dem.df, aes(x, y, fill = elev)) +
  scale_fill_viridis_c("Elevation (m)", na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"))


setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")








###################################
### Merge covars as RasterStack ###
###################################

covars<- stack(dem, green, wet)




#######################################################
### Extract values from raster layer for each track ###
#######################################################
plan(multisession)
path<- extract.covars(data = dat, layers = covars, state.col = "z.post.thresh",
                      dyn_names = c("green","wet"), ind = "season")
#takes 3 min to run

future:::ClusterRegistry("stop")  #close all threads and memory used




#############################################
### Explore relationships among variables ###
#############################################

PerformanceAnalytics::chart.Correlation(path[,2:4])  #no strong corrs (< 0.6)


###################
### Export data ###
###################

# write.csv(path, "Giant Armadillo Resistance Data.csv", row.names = F)


