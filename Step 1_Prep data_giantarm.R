### Load armadillo telemetry data

library(bayesmove)
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


# Calculate step length (also provides turning angle, NSD, and dt)
dat<- prep_data(dat = dat, coord.names = c('x','y'), id = "id")



#######################################
### Import Environmental Covariates ###
#######################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")


## EVI

evi<- brick('GiantArm_evi_monthly.grd')
evi<- crop(evi, extent(dat %>% 
                         summarize(xmin = min(x) - 3000,
                                   xmax = max(x) + 3000,
                                   ymin = min(y) - 3000,
                                   ymax = max(y) + 3000) %>% 
                         unlist()))
evi[getValues(evi) > 1 | getValues(evi) < -1]<- NA  #mask pixels where values are outside of accepted range
evi.s<- scale(evi)


# View distribs of EVI by month
evi.df<- as.data.frame(evi, xy = T)
evi.df<- pivot_longer(evi.df, cols = -c(x,y), names_to = "month", values_to = "evi") %>% 
  mutate_at("month", factor, levels = names(evi))

#density plot
ggplot(evi.df, aes(evi, fill = month)) +
  geom_density(alpha = 0.5) +
  xlim(0,1) +
  labs(x = "EVI", y = "Density") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

#boxplot
ggplot(evi.df, aes(evi, month, color = month)) +
  geom_boxplot(alpha = 0.5) +
  # xlim(0,1) +
  labs(x = "EVI", y = "Month") +
  scale_color_viridis_d(guide = F) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))

#map
ggplot(evi.df, aes(x, y, fill = evi)) +
  geom_tile() +
  labs(x = "Easting", y = "Northing") +
  scale_fill_gradientn("EVI",
                       colors = c('#FFFFFF','#CE7E45','#DF923D','#F1B555','#FCD163','#99B718',
                                  '#74A901','#66A000','#529400','#3E8601','#207401','#056201',
                                  '#004C00','#023B01','#012E01','#011D01','#011301'),
                       na.value = "transparent", limits = c(0,1)) +
  coord_equal() +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ month)








###################################
### Merge covars as RasterStack ###
###################################

# covars<- stack(dem, green, wet)




#######################################################
### Extract values from raster layer for each track ###
#######################################################
plan(multisession)
path<- extract.covars(data = dat, layers = evi, state.col = "z.post.thresh",
                      dyn_names = 'evi', ind = "month")
#takes 2.2 min to run

future:::ClusterRegistry("stop")  #close all threads and memory used




#############################################
### Explore relationships among variables ###
#############################################

# PerformanceAnalytics::chart.Correlation(path[,2:4])  #no strong corrs (< 0.6)


###################
### Export data ###
###################

# write.csv(path, "Giant Armadillo Resistance Data.csv", row.names = F)


