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
dat<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)
dat$date<- as_datetime(dat$date)

# Filter out observations where coords are NA
dat<- dat %>% 
  filter(!is.na(x))


dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))


#########################
### Import LU/LC data ###
#########################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

tmp<- raster('cheia_UTM1.tif')
tmp1<- raster('cheiann_UTM1.tif')


### Modify LU/LC TIFFs 

#Viz rasters
plot(tmp)
plot(tmp1)


#Rename LU/LC rasters
lulc<- tmp1
lulc<- asFactor(lulc)




##########################
### Import NDVI layers ###
##########################


#load files
files<- list.files(getwd(), pattern = "*.grd$")
ndvi.filenames<- files[grep("ndvi", files)]
# ndvi.month<- brick(ndvi.filenames[1])
ndvi.season<- brick(ndvi.filenames[2])

#change extent and dimensions of RasterBricks using resample()
ndvi.season<- resample(ndvi.season, lulc, method = "bilinear")
compareRaster(lulc, ndvi.season)
plot(ndvi.season[[1]]); points(dat$x, dat$y)




##########################
### Import NDWI layers ###
##########################


#load files
files<- list.files(getwd(), pattern = "*.grd$")
ndwi.filenames<- files[grep("ndwi", files)]
# ndwi.month<- brick(ndwi.filenames[1])
ndwi.season<- brick(ndwi.filenames[2])

#change extent and dimensions of RasterBricks using resample()
ndwi.season<- resample(ndwi.season, lulc, method = "bilinear")
compareRaster(lulc, ndwi.season)
plot(ndwi.season[[1]]); points(dat$x, dat$y)





###########################################
### Merge static covars as RasterBricks ###
###########################################

covars<- stack(lulc, ndvi.season, ndwi.season)
names(covars)[1]<- c("lulc")




#######################################################
### Extract values from raster layer for each track ###
#######################################################
# foo<- dat %>% filter(id  %in% c('sara','tex'))
plan(multisession)
path<- extract.covars(data = dat, layers = covars, state.col = NULL, which_cat = "lulc",
                      dyn_names = c("ndvi","ndwi"), ind = "season")
#takes 11 min to run

future:::ClusterRegistry("stop")  #close all threads and memory used




#############################################
### Explore relationships among variables ###
#############################################

PerformanceAnalytics::chart.Correlation(path[,2:6])  #no strong corrs


###################
### Export data ###
###################

# write.csv(path, "Giant Armadillo Resistance Data.csv", row.names = F)


