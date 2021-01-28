
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





## Make spatiotemporal predictions

#extract beta coeffs (mean)
betas<- dat.summ$mean


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
id1<- unique(dat.summ$id)

for (j in 1:nlayers(green)) {
  print(names(green)[j])
  
  tmp<- list()
  
  cov.mat<- cbind(#elev = scale(values(dem)),
                  green = scale(values(green[[j]])),
                  wet = scale(values(wet[[j]])))
  
  for (i in 1:length(id1)) {
    print(i)
    
    resistSurf<- green[[1]]
    raster::values(resistSurf)<- exp(cov.mat %*% dat.summ[which(dat.summ$id == id1[i]), "mean"])
    # resistSurf<- resistSurf * 60  #convert from min to sec
    
    #create as data frame
    resistSurf.df<- as.data.frame(resistSurf, xy=T) %>%
      mutate(id = i)
    names(resistSurf.df)[3]<- "time"
    
    tmp[[i]]<- resistSurf.df
  }
  
  names(tmp)<- names(sort(table(path$id), decreasing = TRUE))
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
                       na.value = "transparent", limits = c(0,15)) +
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


tmp<- bayesmove::df_to_list(dat, "id") %>% 
  purrr::map(., . %>% 
               dplyr::select(season) %>% 
               # unique() %>% 
               unlist())

weights<- tmp %>% 
  map(., ~{table(.x)/length(.x)}) %>% 
  map(., ~{.x[.x !=0]})


#calculate weighted mean resistance across seasons (by N per season)
resist.mean.id<- resist.dyn.df %>% 
  bayesmove::df_to_list(., "id") %>% 
  map2(., tmp, ~{.x %>% 
      filter(season %in% unique(.y))}) %>% 
  map(., bayesmove::df_to_list, "season") %>% 
  map_depth(., 2, pluck, "time")

resist.mean.id2<- resist.mean.id %>% 
  map(., bind_cols) %>% 
  map2(., weights, ~apply(.x, 1, function(x) weighted.mean(x, .y, na.rm = TRUE))) %>% 
  unlist()

resist.mean.id3<- cbind(resist.dyn$Fall$gala[,c("x","y")], time = resist.mean.id2) %>% 
  mutate(id = rep(names(resist.mean.id), each = ncell(green)), .before = "x")
# resist.mean.id3$id<- factor(resist.mean.id3$id, levels = names(resist.mean.id))



## Mean across seasons
ggplot() +
  geom_raster(data = resist.mean.id3, aes(x, y, fill = time)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent", limits = c(0,10)) +
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
  facet_wrap(~ id, nrow = 1)

# ggsave("Giant Armadillo Time Resistance_ID_facet.png", width = 8.5, height = 4,
#        units = "in", dpi = 300)




#calculate weighted mean and variance across IDs (weighted by N for each ID)
resist.pop<- resist.mean.id3 %>% 
  bayesmove::df_to_list(., "id") %>% 
  map_depth(., 1, pluck, "time")

weights2<- table(dat$id)/nrow(dat)
weights2<- weights2[c(1,3,5,4,2,7,6)]  #reorder to match column order of resist.pop

resist.pop.mean<- resist.pop %>% 
  bind_cols() %>% 
  apply(., 1, function(x)  weighted.mean(x, weights2, na.rm = TRUE))

resist.pop.var<- resist.pop %>% 
  bind_cols() %>% 
  apply(., 1, var, na.rm = TRUE)

resist.pop2<- cbind(resist.dyn$Fall$gala[,c("x","y")], mu = resist.pop.mean,
                    sig = resist.pop.var)



## Mean across IDs
ggplot() +
  geom_raster(data = resist.pop2, aes(x, y, fill = mu)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent", limits = c(0,10)) +
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

# ggsave("Giant Armadillo Time Resistance_mean.png", width = 9, height = 5,
#        units = "in", dpi = 300)


## Variance across IDs
ggplot() +
  geom_raster(data = resist.pop2, aes(x, y, fill = sig)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "viridis",
                       na.value = "transparent", limits = c(0,50)) +
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

# ggsave("Giant Armadillo Time Resistance_var.png", width = 9, height = 5,
#        units = "in", dpi = 300)



### Export results


# write.csv(resist.pop2, "Giant Armadillo Resistance summary results.csv", row.names = F)
