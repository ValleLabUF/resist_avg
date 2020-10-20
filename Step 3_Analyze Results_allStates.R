### Visualize Results from Best Models ###

library(ggridges)
library(raster)

### North Pantanal

#look at betas (convert to data frame)
store.betas_N<- data.frame(mod_N$betas[(nburn+1):ngibbs, ])
names(store.betas_N)<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
store.betas.long_N<- tidyr::pivot_longer(store.betas_N,
                                               cols = names(store.betas_N),
                                               names_to = "betas")
store.betas.long_N$betas<- factor(store.betas.long_N$betas,
                                   levels = names(store.betas_N))

ggplot(store.betas.long_N, aes(x=betas, y=value)) +
  geom_boxplot(color="firebrick") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "North Pantanal") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(store.betas.long_N, aes(y=betas, x=value, fill = betas)) +
  geom_density_ridges() +
  scale_fill_viridis_d("Coeffs", guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(y="Beta Coefficients", x="Value", title = "North Pantanal") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))








### South Pantanal

#look at betas (convert to data frame)
store.betas_S<- data.frame(mod_S$betas[(nburn+1):ngibbs, ])
names(store.betas_S)<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
store.betas.long_S<- tidyr::pivot_longer(store.betas_S,
                                               cols = names(store.betas_S),
                                               names_to = "betas")
store.betas.long_S$betas<- factor(store.betas.long_S$betas,
                                        levels = names(store.betas_S))

ggplot(store.betas.long_S, aes(x=betas, y=value)) +
  geom_boxplot(color="darkturquoise") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "South Pantanal") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(store.betas.long_S,aes(y=betas, x=value, fill = betas)) +
  geom_density_ridges() +
  scale_fill_viridis_d("Coeffs", guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(y="Beta Coefficients", x="Value", title = "South Pantanal") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))









#### Create Predictive Surfaces (lulc, t.ar, rain) ####

#Load data to plot tracks
dat<- read.csv("Armadillo HMM Results.csv", header = T, sep = ",")

dat<- dat %>% 
  rename(id = ID) %>% 
  mutate_at("id", as.character) %>% 
  mutate_at("state", as.factor) %>% 
  mutate_at("state", ~recode(., '1' = "Burrow",
                             '2' = "Foraging", '3' = "Transit"))
dat$date<- lubridate::as_datetime(dat$date)

# Separate tracks by region (N or S)
dat.N<- dat %>% filter(region == "N")
dat.S<- dat %>% filter(region == "S")

## North Pantanal

#extract beta coeffs (mean)
betas_N<- colMeans(store.betas_N)


#Need to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulcN<- raster('lulc_N.tif')
lulcN<- as.factor(lulcN)


scaled_t.ar_N<- path.N$t.ar %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))

scaled_rain_N<- path.N$rain %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))


## Foraging ##

##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
xmat<- data.matrix(path.N[,ind])

lulcN.mat<- model.matrix(~factor(getValues(lulcN)) + 0)
colnames(lulcN.mat)<- ind[1:6]


#Min recorded temperature
N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$min, rain = scaled_rain_N$mean)

resistSurfN_minTemp<- lulcN
raster::values(resistSurfN_minTemp)<- exp(N.mat %*% betas_N)
resistSurfN_minTemp<- resistSurfN_minTemp * 60  #convert from min to sec
resistSurfN_minTemp.df<- as.data.frame(resistSurfN_minTemp, xy=T) %>% 
  mutate(temp.level = "Min")


#Avg recorded temperature
N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$mean, rain = scaled_rain_N$mean)

resistSurfN_avgTemp<- lulcN
raster::values(resistSurfN_avgTemp)<- exp(N.mat %*% betas_N)
resistSurfN_avgTemp<- resistSurfN_avgTemp * 60  #convert from min to sec
resistSurfN_avgTemp.df<- as.data.frame(resistSurfN_avgTemp, xy=T) %>% 
  mutate(temp.level = "Avg")

#Max recorded temperature
N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$max, rain = scaled_rain_N$mean)

resistSurfN_maxTemp<- lulcN
raster::values(resistSurfN_maxTemp)<- exp(N.mat %*% betas_N)
resistSurfN_maxTemp<- resistSurfN_maxTemp * 60  #convert from min to sec
resistSurfN_maxTemp.df<- as.data.frame(resistSurfN_maxTemp, xy=T) %>% 
  mutate(temp.level = "Max")


#Combine all results together for each level of temperature
resistSurfN.df<- rbind(resistSurfN_minTemp.df, resistSurfN_avgTemp.df,
                              resistSurfN_maxTemp.df)
resistSurfN.df$temp.level<- factor(resistSurfN.df$temp.level,
                                          levels = c("Min", "Avg", "Max"))
names(resistSurfN.df)[3]<- "time"



ggplot() +
  geom_raster(data = resistSurfN.df, aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time Spent\nper Cell (sec)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "North Pantanal Resistance Surface") +
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
  facet_wrap(~ temp.level)










## South Pantanal

#extract beta coeffs (mean)
betas_S<- colMeans(store.betas_S)

#Seed to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulcS<- raster('lulc_S.tif')
lulcS<- as.factor(lulcS)



scaled_t.ar_S<- path.S$t.ar %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))

scaled_rain_S<- path.S$rain %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))


##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S[,ind])

lulcS.mat<- model.matrix(~factor(getValues(lulcS)) + 0)
colnames(lulcS.mat)<- ind[1:5]


#Min recorded temperature
S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$min, rain = scaled_rain_S$mean)

resistSurfS_minTemp<- lulcS
raster::values(resistSurfS_minTemp)<- exp(S.mat %*% betas_S)
resistSurfS_minTemp<- resistSurfS_minTemp * 60  #convert from min to sec
resistSurfS_minTemp.df<- as.data.frame(resistSurfS_minTemp, xy=T) %>% 
  mutate(temp.level = "Min")

#replace values for water class with NA
resistSurfS_minTemp.df[which(values(lulcS) == 4), 3]<- NA



#Avg recorded temperature
S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$mean, rain = scaled_rain_S$mean)

resistSurfS_avgTemp<- lulcS
raster::values(resistSurfS_avgTemp)<- exp(S.mat %*% betas_S)
resistSurfS_avgTemp<- resistSurfS_avgTemp * 60  #convert from min to sec
resistSurfS_avgTemp.df<- as.data.frame(resistSurfS_avgTemp, xy=T) %>% 
  mutate(temp.level = "Avg")

#replace values for water class with NA
resistSurfS_avgTemp.df[which(values(lulcS) == 4), 3]<- NA



#Max recorded temperature
S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$max, rain = scaled_rain_S$mean)

resistSurfS_maxTemp<- lulcS
raster::values(resistSurfS_maxTemp)<- exp(S.mat %*% betas_S)
resistSurfS_maxTemp<- resistSurfS_maxTemp * 60  #convert from min to sec
resistSurfS_maxTemp.df<- as.data.frame(resistSurfS_maxTemp, xy=T) %>% 
  mutate(temp.level = "Max")

#replace values for water class with NA
resistSurfS_maxTemp.df[which(values(lulcS) == 4), 3]<- NA



#Combine all results together for each level of temperature
resistSurfS.df<- rbind(resistSurfS_minTemp.df, resistSurfS_avgTemp.df,
                              resistSurfS_maxTemp.df)
resistSurfS.df$temp.level<- factor(resistSurfS.df$temp.level,
                                          levels = c("Min", "Avg", "Max"))
names(resistSurfS.df)[3]<- "time"



ggplot() +
  geom_raster(data = resistSurfS.df, aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time Spent\nper Cell (sec)", option = "inferno",
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
  facet_wrap(~ temp.level)

