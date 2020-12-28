
### Visualize Model Results ###

library(tidyverse)
library(lubridate)
library(raster)




#### Create Predictive Surfaces (lulc, ndvi) ####
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




# Load resistance data
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

path<- read.csv('Giant Armadillo Resistance Data.csv', as.is=T)
path$dt<- path$dt/60  #convert to min from sec
path$date<- as_datetime(path$date)


# Filter data for only steps with 5 >= dt >= 9 min
path<- path[path$dt >= 5 & path$dt <= 9 & !is.na(path$dt),]


# Remove state col and rows where is.na(NDVI)
path<- path %>% 
  dplyr::select(-state) %>% 
  filter(!is.na(ndvi))


# Center and scale covariates 
path$ndvi<- as.numeric(scale(path$ndvi, center = TRUE, scale = TRUE))





#Load model results
dat.summ<- read.csv("Giant Armadillo Resistance Results.csv", as.is = T)





## Make spatiotemporal predictions

#extract beta coeffs (mean)
betas<- dat.summ$mean


#Need to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulc<- raster('cheiann_UTM1.tif')
# lulc<- as.factor(lulc)

ndvi.filenames<- list.files(getwd(), pattern = "*.grd$")
ndvi<- brick(ndvi.filenames[2])

#change extent and dimensions of RasterBricks using resample()
ndvi<- resample(ndvi, lulc, method = "bilinear")
compareRaster(lulc, ndvi)



##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("ndvi","X1", "X2", "X3", "X4")
xmat<- data.matrix(path[,ind])

x<- factor(getValues(lulc))
lulc.mat<- model.matrix.lm(~x + 0, na.action = "na.pass")
colnames(lulc.mat)<- ind[2:5]





### Calculate Resistance Surface by ID ###

resist.dyn<- list()
id1<- unique(dat.summ$id)

for (j in 1:nlayers(ndvi)) {
  print(names(ndvi)[j])
  
  tmp<- list()
  
  # Using avg temperature and rainfall estimates
  cov.mat<- cbind(ndvi = values(ndvi[[j]]), lulc.mat)
  
  for (i in 1:length(id1)) {
    print(i)
    
    resistSurf<- lulc
    raster::values(resistSurf)<- exp(cov.mat %*% dat.summ[which(dat.summ$id == id1[i]), "mean"])
    # resistSurf<- resistSurf * 60  #convert from min to sec
    
    #replace values for water class with NA
    # values(resistSurf)[which(values(lulc) == 4)]<- NA
    
    #create as data frame
    resistSurf.df<- as.data.frame(resistSurf, xy=T) %>%
      mutate(id = i)
    names(resistSurf.df)[3]<- "time"
    
    tmp[[i]]<- resistSurf.df
  }
  
  names(tmp)<- names(sort(table(path$id), decreasing = TRUE))
  resist.dyn[[j]]<- tmp
  names(resist.dyn)[j]<- names(ndvi)[j]
}



#Combine all results together for each level of temperature
# resistSurfN.df<- rbind(resistSurfN_minTemp.df, resistSurfN_avgTemp.df,
#                               resistSurfN_maxTemp.df)
# resistSurfN.df$temp.level<- factor(resistSurfN.df$temp.level,
#                                           levels = c("Min", "Avg", "Max"))

resist.dyn.df<- map(resist.dyn, ~bind_rows(., .id = "id")) %>% 
  bind_rows(., .id = "season")


## By ID and season
ggplot() +
  geom_raster(data = resist.dyn.df, aes(x, y, fill = time)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Giant Armadillo Resistance Surface") +
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
  facet_wrap(season ~ id, nrow = 4)

# ggsave("Giant Armadillo Time Resistance_IDxSeason_facet.png", width = 8.5, height = 8, units = "in", dpi = 300)


#calculate mean resistance across IDs
resist.mean.season<- resist.dyn.df %>% 
  bayesmove::df_to_list(., "season") %>% 
  map(., bayesmove::df_to_list, "id") %>% 
  map_depth(., 2, pluck, "time")

resist.mean.season2<- resist.mean.season %>% 
  map(., bind_cols) %>% 
  map(., rowMeans) %>% 
  unlist()

resist.mean.season3<- cbind(resist.dyn$Fall$gala[,c("x","y")], time = resist.mean.season2) %>% 
  mutate(season = rep(c("Fall","Winter","Spring","Summer"), each = ncell(ndvi)), .before = "x")



## Mean across IDs
ggplot() +
  geom_raster(data = resist.mean.season3, aes(x, y, fill = time)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Giant Armadillo Resistance Surface") +
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
  facet_wrap(~ season, nrow = 1)

# ggsave("N Pantanal Time Resistance_mean.png", width = 8.5, height = 8, units = "in", dpi = 300)


## LU/LC map for reference
lulcN.df<- resist.mean.N.df %>% 
  rename(lulc = time) %>% 
  mutate(lulc = values(lulcN))

ggplot() +
  geom_raster(data = lulcN.df, aes(x, y, fill = factor(lulc))) +
  scale_fill_brewer("", palette = "Accent", na.value = "transparent",
                    labels = c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "North Pantanal Land Cover") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))









## South Pantanal

#extract beta coeffs (mean)
# betas_S<- unlist(MCMCpstr(res.S, params = "betas"))
betas_S<- dat.S.summ$mean


#Seed to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulcS<- raster('lulc_S.tif')
# lulcS<- as.factor(lulcS)
ndviS<- raster('ndvi_S.tif')
ndviS<- scale(ndviS, center = T, scale = T)


scaled_t.ar_S<- path.S$t.ar %>% 
  # scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))

scaled_rain_S<- path.S$rain %>% 
  # scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))





##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("ndvi","Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S[,ind])

lulcS.mat<- model.matrix(~factor(getValues(lulcS)) + 0)
colnames(lulcS.mat)<- ind[2:6]


# #Min recorded temperature
# S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$min, rain = scaled_rain_S$mean)
# 
# resistSurfS_minTemp<- lulcS
# raster::values(resistSurfS_minTemp)<- exp(S.mat %*% betas_S)
# resistSurfS_minTemp<- resistSurfS_minTemp * 60  #convert from min to sec
# resistSurfS_minTemp.df<- as.data.frame(resistSurfS_minTemp, xy=T) %>% 
#   mutate(temp.level = "Min")
# 
# #replace values for water class with NA
# resistSurfS_minTemp.df[which(values(lulcS) == 4), 3]<- NA
# 
# 
# 
# #Avg recorded temperature
# S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$mean, rain = scaled_rain_S$mean)
# 
# resistSurfS_avgTemp<- lulcS
# raster::values(resistSurfS_avgTemp)<- exp(S.mat %*% betas_S)
# resistSurfS_avgTemp<- resistSurfS_avgTemp * 60  #convert from min to sec
# resistSurfS_avgTemp.df<- as.data.frame(resistSurfS_avgTemp, xy=T) %>% 
#   mutate(temp.level = "Avg")
# 
# #replace values for water class with NA
# resistSurfS_avgTemp.df[which(values(lulcS) == 4), 3]<- NA
# 
# 
# 
# #Max recorded temperature
# S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$max, rain = scaled_rain_S$mean)
# 
# resistSurfS_maxTemp<- lulcS
# raster::values(resistSurfS_maxTemp)<- exp(S.mat %*% betas_S)
# resistSurfS_maxTemp<- resistSurfS_maxTemp * 60  #convert from min to sec
# resistSurfS_maxTemp.df<- as.data.frame(resistSurfS_maxTemp, xy=T) %>% 
#   mutate(temp.level = "Max")
# 
# #replace values for water class with NA
# resistSurfS_maxTemp.df[which(values(lulcS) == 4), 3]<- NA




### Calculate Resistance Surface by ID ###

resist.S<- list()

# Using avg temperature and rainfall estimates
S.mat<- cbind(ndvi = values(ndviS), lulcS.mat, t.ar = scaled_t.ar_S$mean,
              rain = scaled_rain_S$mean)

for (i in 2:11) {
  print(i)
  
  resistSurfS<- lulcS
  raster::values(resistSurfS)<- exp(S.mat %*% betas_S[12:19] + betas_S[i])
  # resistSurfS<- resistSurfS * 60  #convert from min to sec
  
  #replace values for water class with NA
  values(resistSurfS)[which(values(lulcS) == 4)]<- NA
  
  #create as data frame
  resistSurfS.df<- as.data.frame(resistSurfS, xy=T) %>%
    mutate(id = i)
  
  
  resist.S[[i-1]]<- resistSurfS.df
}



#Combine all results together for each level of temperature
# resistSurfS.df<- rbind(resistSurfS_minTemp.df, resistSurfS_avgTemp.df,
#                               resistSurfS_maxTemp.df)
# resistSurfS.df$temp.level<- factor(resistSurfS.df$temp.level,
#                                           levels = c("Min", "Avg", "Max"))

names(resist.S)<- names(sort(table(path.S$id), decreasing = TRUE))
resist.S.df<- bind_rows(resist.S, .id = "id")
names(resist.S.df)[3]<- "time"


ggplot() +
  geom_raster(data = resist.S.df,
              aes(x, y, fill = time)) +
  geom_path(data = dat.S,
            aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Resistance Surface") +
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
  facet_wrap(~ id)

# ggsave("S Pantanal Time Resistance_facet.png", width = 9, height = 5, units = "in", dpi = 300)



#calculate mean resistance across IDs
resist.mean.S<- resist.S.df %>% 
  group_by(id) %>% 
  group_split() %>% 
  map(., . %>% dplyr::select(time)) %>% 
  bind_cols() %>% 
  rowMeans()
resist.mean.S.df<- cbind(resist.S[[1]][,c("x","y")], time = resist.mean.S)



## Mean across IDs
ggplot() +
  geom_raster(data = resist.mean.S.df, aes(x, y, fill = time)) +
  geom_path(data = dat.S, aes(x, y, group = id), alpha = 0.5, color = "blue") +
  scale_fill_viridis_c("Time Spent\nper Cell (min)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Resistance Surface") +
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

# ggsave("S Pantanal Time Resistance_mean.png", width = 9, height = 5, units = "in", dpi = 300)


## LU/LC map for reference
lulcS.df<- resist.mean.S.df %>% 
  rename(lulc = time) %>% 
  mutate(lulc = values(lulcS))

ggplot() +
  geom_raster(data = lulcS.df, aes(x, y, fill = factor(lulc))) +
  scale_fill_brewer("", palette = "Accent", na.value = "transparent",
                    labels = c("Field", "Forest", "Water", "Pasture", "Road")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Land Cover") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

