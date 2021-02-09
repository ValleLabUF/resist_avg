
### Visualize Model Results ###

library(tidyverse)
library(lubridate)
library(raster)




#### Create Predictive Surfaces (greenness, wetness) ####
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




# Load resistance data
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

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
  mutate(across(elev:wet, scale))




#Load model results
dat.summ<- read.csv("Giant Armadillo Resistance Results.csv", as.is = T)
dat.summ2<- dat.summ[order(dat.summ$id),]




## Make spatiotemporal predictions


#Need to center and scale raster values so comparable to beta coeffs

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")

#Tasseled Cap Greenness
green<- brick('GiantArm_tcgreen_season.grd')
green.df<- as.data.frame(green, xy = T)
green.df2<- pivot_longer(green.df, cols = -c(x,y), names_to = "season", values_to = "green")
green.df2$season<- factor(green.df2$season, levels = names(green))


#Tasseled Cap Wetness
wet<- brick('GiantArm_tcwet_season.grd')
compareRaster(green, wet)

wet.df<- as.data.frame(wet, xy = T)
wet.df2<- pivot_longer(wet.df, cols = -c(x,y), names_to = "season", values_to = "wet")
wet.df2$season<- factor(wet.df2$season, levels = names(wet))


#Elevation
# dem<- raster('giantarm_dem.tif')
# names(dem)<- 'elev'
# dem<- resample(dem, green, method = "bilinear")
# compareRaster(green, dem)
# dem.df<- as.data.frame(dem, xy = TRUE)



setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")


##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("green", "wet")
xmat<- data.matrix(path[,ind])

# x<- factor(getValues(lulc))
# lulc.mat<- model.matrix.lm(~x + 0, na.action = "na.pass")
# colnames(lulc.mat)<- ind[2:5]





### Calculate Resistance Surface by ID ###

resist.dyn<- list()
id1<- unique(dat.summ2$id)

for (j in 1:nlayers(green)) {
  print(names(green)[j])
  
  tmp<- list()
  
  cov.mat<- cbind(1,
                  green = scale(values(green[[j]])),
                  wet = scale(values(wet[[j]])))
  
  for (i in 1:length(id1)) {
    print(i)
    
    resistSurf<- green[[1]]
    raster::values(resistSurf)<- exp(cov.mat %*% dat.summ2[which(dat.summ2$id == id1[i]), "mean"])
    # resistSurf<- resistSurf * 60  #convert from min to sec
    
    #create as data frame
    resistSurf.df<- as.data.frame(resistSurf, xy=T) %>%
      mutate(id = i)
    names(resistSurf.df)[3]<- "time"
    
    tmp[[i]]<- resistSurf.df
  }
  
  # names(tmp)<- names(sort(table(path$id), decreasing = TRUE))
  names(tmp)<- id1
  resist.dyn[[j]]<- tmp
  names(resist.dyn)[j]<- names(green)[j]
}



#Combine all results together for each level of temperature
# resistSurfN.df<- rbind(resistSurfN_minTemp.df, resistSurfN_avgTemp.df,
#                               resistSurfN_maxTemp.df)
# resistSurfN.df$temp.level<- factor(resistSurfN.df$temp.level,
#                                           levels = c("Min", "Avg", "Max"))

resist.dyn.df<- map(resist.dyn, ~bind_rows(., .id = "id")) %>% 
  bind_rows(., .id = "season")
resist.dyn.df$season<- factor(resist.dyn.df$season, levels = names(green))


## By ID and season
ggplot() +
  geom_raster(data = resist.dyn.df, aes(x, y, fill = time)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent", limits = c(0,7)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Resistance Surface") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1)) +
  facet_grid(season ~ id)

# ggsave("Giant Armadillo Time Resistance_IDxSeason_facet.png", width = 9, height = 7.5,
#        units = "in", dpi = 300)


#calculate mean resistance at population level (using mean (of means) for ID coeff estimates)
resist.mean<- list()
mean1<- dat.summ2 %>% 
  group_by(coeff) %>% 
  summarize(mean = mean(mean)) %>% 
  ungroup() %>% 
  dplyr::select(mean) %>% 
  unlist()
mean1<- mean1[c(2,1,3)]  #reorder to match cov.mat (int, green, wet)

for (j in 1:nlayers(green)) {
  print(names(green)[j])
  
  cov.mat<- cbind(1,
                  green = scale(values(green[[j]])),
                  wet = scale(values(wet[[j]])))
  
  
  resistSurf<- green[[1]]
  raster::values(resistSurf)<- exp(cov.mat %*% mean1)
  
  #create as data frame
  resistSurf.df<- as.data.frame(resistSurf, xy=T)
  names(resistSurf.df)[3]<- "time"
  
  resist.mean[[j]]<- resistSurf.df
  names(resist.mean)[j]<- names(green)[j]
  
}


resist.mean.df<- bind_rows(resist.mean, .id = "season")
resist.mean.df$season<- factor(resist.mean.df$season, levels = names(green))



## Mean across IDs
ggplot() +
  geom_raster(data = resist.mean.df, aes(x, y, fill = time)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Resistance Surface") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1)) +
  facet_wrap(~ season, nrow = 2)

# ggsave("Giant Armadillo Time Resistance_ID_facet.png", width = 8.5, height = 4,
#        units = "in", dpi = 300)



# Fall ONLY data (to send to Denis)
resist.mean.fall<- resist.mean.df %>% 
  filter(season == "Fall")




### Export results


# write.csv(resist.mean.df, "Giant Armadillo Resistance summary results.csv", row.names = F)
# write.csv(resist.mean.fall, "Giant Armadillo Resistance summary results_Fall.csv", row.names = F)
