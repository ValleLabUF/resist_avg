### Visualize Results from Best Models ###

library(ggridges)
library(raster)

### North Pantanal

## Foraging

# #look at cross-correlation among betas
# k=cor(cbind(store.b.forage_N,store.betas.forage_N))
# k[k < 0.5 & k > -0.5]=NA
# diag(k)=NA
# k

#look at betas (convert to data frame)
store.betas_Nforage<- data.frame(mod.forage_N$betas[(nburn+1):ngibbs, ])
names(store.betas_Nforage)<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
store.betas.long_Nforage<- tidyr::pivot_longer(store.betas_Nforage,
                                               cols = names(store.betas_Nforage),
                                               names_to = "betas")
store.betas.long_Nforage$betas<- factor(store.betas.long_Nforage$betas,
                                   levels = names(store.betas_Nforage))

ggplot(store.betas.long_Nforage, aes(x=betas, y=value)) +
  geom_boxplot(color="firebrick") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "North Pantanal (Foraging)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# # w/o intercept
# ggplot(store.betas.long_Nforage %>% filter(betas != "int"), aes(x=betas, y=value)) +
#   geom_boxplot(color="firebrick") +
#   geom_hline(yintercept = 0, size = 0.5) +
#   labs(x="Beta Coefficients", y="Value", title = "North Pantanal (Foraging)") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.title = element_text(size = 12))

ggplot(store.betas.long_Nforage, aes(y=betas, x=value, fill = betas)) +
  geom_density_ridges() +
  scale_fill_viridis_d("Coeffs", guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(y="Beta Coefficients", x="Value", title = "North Pantanal (Foraging)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))




## Transit

#look at betas (convert to data frame)
store.betas_Ntransit<- data.frame(mod.transit_N$betas[(nburn+1):ngibbs, ])
names(store.betas_Ntransit)<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
store.betas.long_Ntransit<- tidyr::pivot_longer(store.betas_Ntransit,
                                               cols = names(store.betas_Ntransit),
                                               names_to = "betas")
store.betas.long_Ntransit$betas<- factor(store.betas.long_Ntransit$betas,
                                        levels = names(store.betas_Ntransit))

ggplot(store.betas.long_Ntransit, aes(x=betas, y=value)) +
  geom_boxplot(color="firebrick") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "North Pantanal (Transit)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# w/o intercept
ggplot(store.betas.long_Ntransit %>% filter(betas != "int"), aes(x=betas, y=value)) +
  geom_boxplot(color="firebrick") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "North Pantanal (Transit)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(store.betas.long_Ntransit %>% filter(betas != "int"), aes(y=betas, x=value,
                                                                 fill = betas)) +
  geom_density_ridges() +
  scale_fill_viridis_d("Coeffs", guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(y="Beta Coefficients", x="Value", title = "North Pantanal (Transit)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))






### South Pantanal

## Foraging

#look at betas (convert to data frame)
store.betas_Sforage<- data.frame(mod.forage_S$betas[(nburn+1):ngibbs, ])
names(store.betas_Sforage)<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
store.betas.long_Sforage<- tidyr::pivot_longer(store.betas_Sforage,
                                               cols = names(store.betas_Sforage),
                                               names_to = "betas")
store.betas.long_Sforage$betas<- factor(store.betas.long_Sforage$betas,
                                        levels = names(store.betas_Sforage))

ggplot(store.betas.long_Sforage, aes(x=betas, y=value)) +
  geom_boxplot(color="darkturquoise") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "South Pantanal (Foraging)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# w/o intercept
ggplot(store.betas.long_Sforage %>% filter(betas != "int"), aes(x=betas, y=value)) +
  geom_boxplot(color="darkturquoise") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "South Pantanal (Foraging)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(store.betas.long_Sforage %>% filter(betas != "int"),aes(y=betas, x=value, fill = betas)) +
  geom_density_ridges() +
  scale_fill_viridis_d("Coeffs", guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(y="Beta Coefficients", x="Value", title = "South Pantanal (Foraging)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))




## Transit

#look at betas (convert to data frame)
store.betas_Stransit<- data.frame(mod.transit_S$betas[(nburn+1):ngibbs, ])
names(store.betas_Stransit)<-  c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
store.betas.long_Stransit<- tidyr::pivot_longer(store.betas_Stransit,
                                                cols = names(store.betas_Stransit),
                                                names_to = "betas")
store.betas.long_Stransit$betas<- factor(store.betas.long_Stransit$betas,
                                         levels = names(store.betas_Stransit))

ggplot(store.betas.long_Stransit, aes(x=betas, y=value)) +
  geom_boxplot(color="darkturquoise") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "South Pantanal (Transit)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# w/o intercept
ggplot(store.betas.long_Stransit %>% filter(betas != "int"), aes(x=betas, y=value)) +
  geom_boxplot(color="darkturquoise") +
  geom_hline(yintercept = 0, size = 0.5) +
  labs(x="Beta Coefficients", y="Value", title = "South Pantanal (Transit)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggplot(store.betas.long_Stransit %>% filter(betas != "int"), aes(y=betas, x=value,
                                                                 fill = betas)) +
  geom_density_ridges() +
  scale_fill_viridis_d("Coeffs", guide = guide_legend(reverse = TRUE)) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(y="Beta Coefficients", x="Value", title = "South Pantanal (Transit)") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12))









#### Create Predictive Surfaces (dist2rd, slope, ndvi, t.ar, rain) ####

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
betas.forage_N<- colMeans(store.betas_Nforage)


#Need to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulcN<- raster('lulc_N.tif')
lulcN<- as.factor(lulcN)


# Mask all unused pixels
ind_N<- unique(cellFromXY(lulcN, dat.N[, c("x","y")]))
covars.N_masked<- lulcN
covars.N_masked[setdiff(1:ncell(covars.N_masked), ind_N)] <- NA


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
xmat<- data.matrix(path.N.forage[,ind])

lulcN.mat<- model.matrix(~factor(getValues(lulcN)) + 0)
colnames(lulcN.mat)<- ind[1:6]


#Min recorded temperature
forageN.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$min, rain = scaled_rain_N$mean)

resistSurfN.forage_minTemp<- lulcN
raster::values(resistSurfN.forage_minTemp)<- exp(forageN.mat %*% betas.forage_N)
resistSurfN.forage_minTemp<- resistSurfN.forage_minTemp * 60  #convert from min to sec
resistSurfN.forage_minTemp.df<- as.data.frame(resistSurfN.forage_minTemp, xy=T) %>% 
  mutate(temp.level = "Min")


#Avg recorded temperature
forageN.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$mean, rain = scaled_rain_N$mean)

resistSurfN.forage_avgTemp<- lulcN
raster::values(resistSurfN.forage_avgTemp)<- exp(forageN.mat %*% betas.forage_N)
resistSurfN.forage_avgTemp<- resistSurfN.forage_avgTemp * 60  #convert from min to sec
resistSurfN.forage_avgTemp.df<- as.data.frame(resistSurfN.forage_avgTemp, xy=T) %>% 
  mutate(temp.level = "Avg")

#Max recorded temperature
forageN.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$max, rain = scaled_rain_N$mean)

resistSurfN.forage_maxTemp<- lulcN
raster::values(resistSurfN.forage_maxTemp)<- exp(forageN.mat %*% betas.forage_N)
resistSurfN.forage_maxTemp<- resistSurfN.forage_maxTemp * 60  #convert from min to sec
resistSurfN.forage_maxTemp.df<- as.data.frame(resistSurfN.forage_maxTemp, xy=T) %>% 
  mutate(temp.level = "Max")


#Combine all results together for each level of lunar illumination
resistSurfN.forage.df<- rbind(resistSurfN.forage_minTemp.df, resistSurfN.forage_avgTemp.df,
                              resistSurfN.forage_maxTemp.df)
resistSurfN.forage.df$temp.level<- factor(resistSurfN.forage.df$temp.level,
                                          levels = c("Min", "Avg", "Max"))
names(resistSurfN.forage.df)[3]<- "time"



forageN_plot<- ggplot() +
  geom_raster(data = resistSurfN.forage.df, aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time Spent\nper Cell (sec)", option = "inferno",
                       na.value = "transparent") +
  # geom_point(data = dat.N %>% filter(state == "Foraging"), aes(x, y, color = id),
  #            size = 0.5, alpha = 0.2, show.legend = F) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "North Pantanal Foraging Resistance Surface") +
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





## Transit ##

##Perform raster math using beta coeffs

#extract beta coeffs (mean)
betas.transit_N<- colMeans(store.betas_Ntransit)

#Min recorded temperature
transitN.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$min, rain = scaled_rain_N$mean)

resistSurfN.transit_minTemp<- lulcN
raster::values(resistSurfN.transit_minTemp)<- exp(transitN.mat %*% betas.transit_N)
resistSurfN.transit_minTemp<- resistSurfN.transit_minTemp * 60  #convert from min to sec
resistSurfN.transit_minTemp.df<- as.data.frame(resistSurfN.transit_minTemp, xy=T) %>% 
  mutate(temp.level = "Min")

#replace values for water class with NA
resistSurfN.transit_minTemp.df[which(values(lulcN) == 4), 3]<- NA


#Avg recorded temperature
transitN.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$mean, rain = scaled_rain_N$mean)

resistSurfN.transit_avgTemp<- lulcN
raster::values(resistSurfN.transit_avgTemp)<- exp(transitN.mat %*% betas.transit_N)
resistSurfN.transit_avgTemp<- resistSurfN.transit_avgTemp * 60  #convert from min to sec
resistSurfN.transit_avgTemp.df<- as.data.frame(resistSurfN.transit_avgTemp, xy=T) %>% 
  mutate(temp.level = "Avg")

#replace values for water class with NA
resistSurfN.transit_avgTemp.df[which(values(lulcN) == 4), 3]<- NA


#Max recorded temperature
transitN.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$max, rain = scaled_rain_N$mean)

resistSurfN.transit_maxTemp<- lulcN
raster::values(resistSurfN.transit_maxTemp)<- exp(transitN.mat %*% betas.transit_N)
resistSurfN.transit_maxTemp<- resistSurfN.transit_maxTemp * 60  #convert from min to sec
resistSurfN.transit_maxTemp.df<- as.data.frame(resistSurfN.transit_maxTemp, xy=T) %>% 
  mutate(temp.level = "Max")

#replace values for water class with NA
resistSurfN.transit_maxTemp.df[which(values(lulcN) == 4), 3]<- NA



#Combine all results together for each level of lunar illumination
resistSurfN.transit.df<- rbind(resistSurfN.transit_minTemp.df, resistSurfN.transit_avgTemp.df,
                              resistSurfN.transit_maxTemp.df)
resistSurfN.transit.df$temp.level<- factor(resistSurfN.transit.df$temp.level,
                                          levels = c("Min", "Avg", "Max"))
names(resistSurfN.transit.df)[3]<- "time"



ggplot() +
  geom_raster(data = resistSurfN.transit.df, aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time Spent\nper Cell (sec)", option = "inferno",
                       na.value = "transparent") +
  # geom_point(data = dat.N %>% filter(state == "Foraging"), aes(x, y, color = id),
  #            size = 0.5, alpha = 0.2, show.legend = F) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "North Pantanal Transit Resistance Surface") +
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
betas.forage_S<- colMeans(store.betas_Sforage)

#Seed to center and scale raster values so comparable to beta coeffs

#Load env raster data
lulcS<- raster('lulc_S.tif')
lulcS<- as.factor(lulcS)


# Mask all unused pixels
ind_S<- unique(cellFromXY(lulcS, dat.S[, c("x","y")]))
covars.S_masked<- lulcS
covars.S_masked[setdiff(1:ncell(covars.S_masked), ind_S)] <- NA


scaled_t.ar_S<- path.S$t.ar %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))

scaled_rain_S<- path.S$rain %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame() %>% 
  summarise(min=min(V1, na.rm = T), mean=mean(V1, na.rm = T), max=max(V1, na.rm = T))


## Foraging ##

##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S.forage[,ind])

lulcS.mat<- model.matrix(~factor(getValues(lulcS)) + 0)
colnames(lulcS.mat)<- ind[1:5]


#Min recorded temperature
forageS.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$min, rain = scaled_rain_S$mean)

resistSurfS.forage_minTemp<- lulcS
raster::values(resistSurfS.forage_minTemp)<- exp(forageS.mat %*% betas.forage_S)
resistSurfS.forage_minTemp<- resistSurfS.forage_minTemp * 60  #convert from min to sec
resistSurfS.forage_minTemp.df<- as.data.frame(resistSurfS.forage_minTemp, xy=T) %>% 
  mutate(temp.level = "Min")


#Avg recorded temperature
forageS.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$mean, rain = scaled_rain_S$mean)

resistSurfS.forage_avgTemp<- lulcS
raster::values(resistSurfS.forage_avgTemp)<- exp(forageS.mat %*% betas.forage_S)
resistSurfS.forage_avgTemp<- resistSurfS.forage_avgTemp * 60  #convert from min to sec
resistSurfS.forage_avgTemp.df<- as.data.frame(resistSurfS.forage_avgTemp, xy=T) %>% 
  mutate(temp.level = "Avg")

#Max recorded temperature
forageS.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$max, rain = scaled_rain_S$mean)

resistSurfS.forage_maxTemp<- lulcS
raster::values(resistSurfS.forage_maxTemp)<- exp(forageS.mat %*% betas.forage_S)
resistSurfS.forage_maxTemp<- resistSurfS.forage_maxTemp * 60  #convert from min to sec
resistSurfS.forage_maxTemp.df<- as.data.frame(resistSurfS.forage_maxTemp, xy=T) %>% 
  mutate(temp.level = "Max")


#Combine all results together for each level of lunar illumination
resistSurfS.forage.df<- rbind(resistSurfS.forage_minTemp.df, resistSurfS.forage_avgTemp.df,
                              resistSurfS.forage_maxTemp.df)
resistSurfS.forage.df$temp.level<- factor(resistSurfS.forage.df$temp.level,
                                          levels = c("Min", "Avg", "Max"))
names(resistSurfS.forage.df)[3]<- "time"



ggplot() +
  geom_raster(data = resistSurfS.forage.df, aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time Spent\nper Cell (sec)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Foraging Resistance Surface") +
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





## Transit ##

##Perform raster math using beta coeffs

#extract beta coeffs (mean)
betas.transit_S<- colMeans(store.betas_Stransit)

#Make predictions using posterior mean of betas
ind<- c("Field", "Forest", "Water", "Pasture", "Road", "t.ar", "rain")
xmat<- data.matrix(path.S.transit[,ind])

lulcS.mat<- model.matrix(~factor(getValues(lulcS)) + 0)
colnames(lulcS.mat)<- ind[1:5]


#Min recorded temperature
transitS.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$min, rain = scaled_rain_S$mean)

resistSurfS.transit_minTemp<- lulcS
raster::values(resistSurfS.transit_minTemp)<- exp(transitS.mat %*% betas.transit_S)
resistSurfS.transit_minTemp<- resistSurfS.transit_minTemp * 60  #convert from min to sec
resistSurfS.transit_minTemp.df<- as.data.frame(resistSurfS.transit_minTemp, xy=T) %>% 
  mutate(temp.level = "Min")


#Avg recorded temperature
transitS.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$mean, rain = scaled_rain_S$mean)

resistSurfS.transit_avgTemp<- lulcS
raster::values(resistSurfS.transit_avgTemp)<- exp(transitS.mat %*% betas.transit_S)
resistSurfS.transit_avgTemp<- resistSurfS.transit_avgTemp * 60  #convert from min to sec
resistSurfS.transit_avgTemp.df<- as.data.frame(resistSurfS.transit_avgTemp, xy=T) %>% 
  mutate(temp.level = "Avg")

#Max recorded temperature
transitS.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$max, rain = scaled_rain_S$mean)

resistSurfS.transit_maxTemp<- lulcS
raster::values(resistSurfS.transit_maxTemp)<- exp(transitS.mat %*% betas.transit_S)
resistSurfS.transit_maxTemp<- resistSurfS.transit_maxTemp * 60  #convert from min to sec
resistSurfS.transit_maxTemp.df<- as.data.frame(resistSurfS.transit_maxTemp, xy=T) %>% 
  mutate(temp.level = "Max")


#Combine all results together for each level of lunar illumination
resistSurfS.transit.df<- rbind(resistSurfS.transit_minTemp.df, resistSurfS.transit_avgTemp.df,
                              resistSurfS.transit_maxTemp.df)
resistSurfS.transit.df$temp.level<- factor(resistSurfS.transit.df$temp.level,
                                          levels = c("Min", "Avg", "Max"))
names(resistSurfS.transit.df)[3]<- "time"



ggplot() +
  geom_raster(data = resistSurfS.transit.df, aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time Spent\nper Cell (sec)", option = "inferno",
                       na.value = "transparent") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "South Pantanal Transit Resistance Surface") +
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
