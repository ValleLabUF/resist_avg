
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


# Remove state col and rows where is.na(NDWI)
path<- path %>% 
  dplyr::select(-state) %>% 
  filter(!is.na(ndwi))


# Center and scale covariates 
path$ndwi<- as.numeric(scale(path$ndwi, center = TRUE, scale = TRUE))





#Load model results
dat.summ<- read.csv("Giant Armadillo Resistance Results.csv", as.is = T)





## Make spatiotemporal predictions

#extract beta coeffs (mean)
betas<- dat.summ$mean


#Need to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulc<- raster('cheiann_UTM1.tif')
# lulc<- as.factor(lulc)

files<- list.files(getwd(), pattern = "*.grd$")
ndwi.filenames<- files[grep("ndwi", files)]
ndwi<- brick(ndwi.filenames[2])

#change extent and dimensions of RasterBricks using resample()
ndwi<- resample(ndwi, lulc, method = "bilinear")
compareRaster(lulc, ndwi)



##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("ndwi","X1", "X2", "X3", "X4")
xmat<- data.matrix(path[,ind])

x<- factor(getValues(lulc))
lulc.mat<- model.matrix.lm(~x + 0, na.action = "na.pass")
colnames(lulc.mat)<- ind[2:5]





### Calculate Resistance Surface by ID ###

resist.dyn<- list()
id1<- unique(dat.summ$id)

for (j in 1:nlayers(ndwi)) {
  print(names(ndwi)[j])
  
  tmp<- list()
  
  cov.mat<- cbind(ndwi = scale(values(ndwi[[j]]), center = T, scale = T), lulc.mat)
  
  for (i in 1:length(id1)) {
    print(i)
    
    resistSurf<- lulc
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
  names(resist.dyn)[j]<- names(ndwi)[j]
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
  facet_grid(season ~ id)

# ggsave("Giant Armadillo Time Resistance_IDxSeason_facet.png", width = 8.5, height = 8,
#        units = "in", dpi = 300)


#calculate mean resistance across seasons
resist.mean.id<- resist.dyn.df %>% 
  bayesmove::df_to_list(., "id") %>% 
  map(., bayesmove::df_to_list, "season") %>% 
  map_depth(., 2, pluck, "time")

resist.mean.id2<- resist.mean.id %>% 
  map(., bind_cols) %>% 
  map(., rowMeans, na.rm = TRUE) %>% 
  unlist()

resist.mean.id3<- cbind(resist.dyn$Fall$gala[,c("x","y")], time = resist.mean.id2) %>% 
  mutate(id = rep(names(resist.mean.id), each = ncell(ndwi)), .before = "x")
# resist.mean.id3$id<- factor(resist.mean.id3$id, levels = names(resist.mean.id))



## Mean across seasons
ggplot() +
  geom_raster(data = resist.mean.id3, aes(x, y, fill = time)) +
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
  facet_wrap(~ id, nrow = 1)

# ggsave("Giant Armadillo Time Resistance_ID_facet.png", width = 8.5, height = 4,
#        units = "in", dpi = 300)




#calculate mean and variance across IDs
resist.pop<- resist.mean.id3 %>% 
  bayesmove::df_to_list(., "id") %>% 
  map_depth(., 1, pluck, "time")

resist.pop.mean<- resist.pop %>% 
  bind_cols() %>% 
  rowMeans(., na.rm = TRUE)

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
                       na.value = "transparent") +
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

# ggsave("Giant Armadillo Time Resistance_var.png", width = 9, height = 5,
#        units = "in", dpi = 300)




## LU/LC map for reference
lulc.df<- resist.mean.id3 %>% 
  rename(lulc = time)
lulc.df$lulc<- raster::values(lulc)

ggplot() +
  geom_raster(data = lulc.df, aes(x, y, fill = factor(lulc))) +
  scale_fill_brewer("", palette = "Accent", na.value = "transparent",
                    labels = c("Forest", "Closed Savanna", "Open Savanna", "Floodable","")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Pantanal Land Cover") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("Giant Armadillo LULC.png", width = 6, height = 5, units = "in", dpi = 300)
