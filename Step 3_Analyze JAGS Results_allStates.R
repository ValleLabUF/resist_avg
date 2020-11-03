### Visualize Model Results ###

library(tidyverse)
library(raster)




#### Create Predictive Surfaces (lulc, t.ar, rain) ####

#Load original data to plot tracks
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





# N and S IDs separated
path.N<- read.csv('N Armadillo Resistance Data_LULC.csv', as.is=T)
path.S<- read.csv('S Armadillo Resistance Data_LULC.csv', as.is=T)
names(path.N)[2:7]<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")
names(path.S)[2:6]<- c("Field", "Forest", "Water", "Pasture", "Road")

path.N$dt<- path.N$dt/60  #convert to min from sec
path.S$dt<- path.S$dt/60

# Filter data for only steps with 3 >= dt >= 7 min
path.N<- path.N[path.N$dt >= 3 & path.N$dt <= 7 & !is.na(path.N$dt),]
path.S<- path.S[path.S$dt >= 3 & path.S$dt <= 7 & !is.na(path.S$dt),]

# Remove Burrow behavior observations
path.N<- path.N %>% 
  filter(state != 'Burrow')
path.S<- path.S %>% 
  filter(state != 'Burrow')

# Center and Scale covariates 
path.N<- path.N %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)

path.S<- path.S %>% 
  mutate_at(c("t.ar","rain"),
            ~scale(., center = TRUE, scale = TRUE)) %>% 
  drop_na(t.ar)





#Load model results
dat.N.summ<- read.csv("N Armadillo Resistance Results.csv")
dat.S.summ<- read.csv("S Armadillo Resistance Results.csv")





## North Pantanal

#extract beta coeffs (mean)
# betas_N<- unlist(MCMCpstr(res.N, params = "betas"))
betas_N<- dat.N.summ$mean


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



##Perform raster math using beta coeffs

#Make predictions using posterior mean of betas
ind<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest", "t.ar", "rain")
xmat<- data.matrix(path.N[,ind])

lulcN.mat<- model.matrix(~factor(getValues(lulcN)) + 0)
colnames(lulcN.mat)<- ind[1:6]


# #Min recorded temperature
# N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$min, rain = scaled_rain_N$mean)
# 
# resistSurfN_minTemp<- lulcN
# raster::values(resistSurfN_minTemp)<- exp(N.mat %*% betas_N)
# resistSurfN_minTemp<- resistSurfN_minTemp * 60  #convert from min to sec
# resistSurfN_minTemp.df<- as.data.frame(resistSurfN_minTemp, xy=T) %>% 
#   mutate(temp.level = "Min")
# 
# #replace values for water class with NA
# resistSurfN_minTemp.df[which(values(lulcN) == 4), 3]<- NA
# 
# 
# 
# #Avg recorded temperature
# N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$mean, rain = scaled_rain_N$mean)
# 
# resistSurfN_avgTemp<- lulcN
# raster::values(resistSurfN_avgTemp)<- exp(N.mat %*% betas_N)
# resistSurfN_avgTemp<- resistSurfN_avgTemp * 60  #convert from min to sec
# resistSurfN_avgTemp.df<- as.data.frame(resistSurfN_avgTemp, xy=T) %>% 
#   mutate(temp.level = "Avg")
# 
# #replace values for water class with NA
# resistSurfN_avgTemp.df[which(values(lulcN) == 4), 3]<- NA
# 
# 
# 
# #Max recorded temperature
# N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$max, rain = scaled_rain_N$mean)
# 
# resistSurfN_maxTemp<- lulcN
# raster::values(resistSurfN_maxTemp)<- exp(N.mat %*% betas_N)
# resistSurfN_maxTemp<- resistSurfN_maxTemp * 60  #convert from min to sec
# resistSurfN_maxTemp.df<- as.data.frame(resistSurfN_maxTemp, xy=T) %>% 
#   mutate(temp.level = "Max")
# 
# #replace values for water class with NA
# resistSurfN_maxTemp.df[which(values(lulcN) == 4), 3]<- NA





### Calculate Resistance Surface by ID ###

resist.N<- list()

# Using avg temperature and rainfall estimates
N.mat<- cbind(lulcN.mat, t.ar = scaled_t.ar_N$mean, rain = scaled_rain_N$mean)

for (i in 2:11) {
  print(i)
  
  resistSurfN<- lulcN
  raster::values(resistSurfN)<- exp(N.mat %*% betas_N[12:19] + betas_N[i])
  resistSurfN<- resistSurfN * 60  #convert from min to sec
  
  #replace values for water class with NA
  values(resistSurfN)[which(values(lulcN) == 4)]<- NA
  
  #create as data frame
  resistSurfN.df<- as.data.frame(resistSurfN, xy=T) %>%
    mutate(id = i)
  
  
  resist.N[[i-1]]<- resistSurfN.df
}


#Combine all results together for each level of temperature
# resistSurfN.df<- rbind(resistSurfN_minTemp.df, resistSurfN_avgTemp.df,
#                               resistSurfN_maxTemp.df)
# resistSurfN.df$temp.level<- factor(resistSurfN.df$temp.level,
#                                           levels = c("Min", "Avg", "Max"))

names(resist.N)<- names(sort(table(path.N$id), decreasing = TRUE))
resist.N.df<- bind_rows(resist.N, .id = "id")
names(resist.N.df)[3]<- "time"


## By ID
ggplot() +
  geom_raster(data = resist.N.df, aes(x, y, fill = time)) +
  geom_path(data = dat.N, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
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
  facet_wrap(~ id)


#calculate mean resistance across IDs
resist.mean.N<- resist.N.df %>% 
  group_by(id) %>% 
  group_split() %>% 
  map(., . %>% dplyr::select(time)) %>% 
  bind_cols() %>% 
  rowMeans()
resist.mean.N.df<- cbind(resist.N[[1]][,c("x","y")], time = resist.mean.N)



## Mean across IDs
ggplot() +
  geom_raster(data = resist.mean.N.df, aes(x, y, fill = time)) +
  geom_path(data = dat.N, aes(x, y, group = id), alpha = 0.5, color = "blue") +
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
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))


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
S.mat<- cbind(lulcS.mat, t.ar = scaled_t.ar_S$mean, rain = scaled_rain_S$mean)

for (i in 2:11) {
  print(i)
  
  resistSurfS<- lulcS
  raster::values(resistSurfS)<- exp(S.mat %*% betas_S[12:18] + betas_S[i])
  resistSurfS<- resistSurfS * 60  #convert from min to sec
  
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

tmp1<- resist.S.df %>% filter(id %in% unique(resist.S.df$id)[1:5])

S.plot1<- ggplot() +
  geom_raster(data = tmp1,
              aes(x, y, fill = time)) +
  geom_path(data = dat.S %>% filter(id %in% unique(resist.S.df$id)[1:5]),
            aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
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
  facet_wrap(~ id)

ggsave("S Pantanal ID-level Resistance Surfaces1.png", plot = S.plot1, width = 9, height = 9,
       units = "in", dpi = 300)



tmp2<- resist.S.df %>% filter(id %in% unique(resist.S.df$id)[6:10])

S.plot2<- ggplot() +
  geom_raster(data = tmp2,
              aes(x, y, fill = time)) +
  geom_path(data = dat.S %>% filter(id %in% unique(resist.S.df$id)[6:10]),
            aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
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
  facet_wrap(~ id)

ggsave("S Pantanal ID-level Resistance Surfaces2.png", plot = S.plot2, width = 9, height = 9,
       units = "in", dpi = 300)



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
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1))


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

