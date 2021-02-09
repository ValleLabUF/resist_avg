### Must first have run Step 3_Analyze Results_giantarm.R and giantarm_ssf.R already

library(tidyverse)
library(lubridate)
library(ggnewscale)
library(raster)

### Load data

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


#Tasseled Cap Greenness
green<- brick('GiantArm_tcgreen_season.grd')

#Tasseled Cap Wetness
wet<- brick('GiantArm_tcwet_season.grd')
compareRaster(green, wet)



setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")
ssf.pop2<- read.csv("Giant Armadillo SSF summary results.csv", as.is = T)
resist.pop2<- read.csv("Giant Armadillo Resistance summary results.csv", as.is = T)





### Plot Velocity-based Resistance vs Habitat Selection

dat.comp<- cbind(ssf.pop2, time = resist.pop2$time)
dat.comp$season<- factor(dat.comp$season, levels = names(green))

quant_99.9<- quantile(na.omit(dat.comp$time), 0.999)
dat.comp2<- dat.comp %>% 
  filter(time <= quant_99.9)
time.mid<- (max(dat.comp2$time, na.rm = T) + min(dat.comp2$time, na.rm = T))/2



ggplot(dat.comp, aes(time, sel, color = season)) +
  geom_point(size = 1, alpha = 0.25) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.mid)) +
  theme_bw() +
  ylim(0,1) +
  xlim(0, max(dat.comp$time, na.rm = T)) +
  labs(x= "Time (min)", y="Habitat Preference") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10)) +
  facet_wrap(~ season, nrow = 2)

# ggsave("Giant Armadillo Quadrant_lulc.png", width = 7, height = 5,
#        units = "in", dpi = 300)



# ggplot(dat.comp2, aes(res, sel, color = wet)) +
#   geom_hline(aes(yintercept = 0.5)) +
#   geom_vline(aes(xintercept = time.mid)) +
#   scale_color_distiller("Wetness", palette = "Blues", direction = 1) +
#   geom_point(na.rm = T) +
#   theme_bw() +
#   ylim(0,1) +
#   labs(x= "Time (min)", y="Habitat Preference") +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 10))

# ggsave("Giant Armadillo Quadrant_ndwi.png", width = 7, height = 5,
#        units = "in", dpi = 300)





### Reclassify pixels and map
dat.comp3<- dat.comp
ind<- which(dat.comp3$time > 7)
dat.comp3[ind, c("sel","time")]<- NA  #set values to NA when time > 7 min
  


## Using 2 continuous color scales
ggplot() +
  geom_raster(data = dat.comp3, aes(x, y, fill = time), alpha = 0.75) +
  scale_fill_gradient("Time Spent\nper Cell (min)", low = "white", high = "royalblue3",
                       na.value = "transparent") +
  new_scale_fill() +
  geom_raster(data = dat.comp3, aes(x, y, fill = sel), alpha = 0.75) +
  scale_fill_gradient("Selection", low = "red4", high = "white",
                       na.value = "transparent", limits = c(0,1)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "black") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))  +
  facet_wrap(~ season)

# ggsave("Giant Armadillo Connectivity_cont.png", width = 7, height = 5,
#        units = "in", dpi = 300)




## Creation of 4 classes

all.dat<- dat.comp3 %>% 
  mutate(fun_class = case_when(.$time > time.mid & .$sel > 0.5 ~ "Slow-Preferred",
                               .$time > time.mid & .$sel < 0.5 ~ "Slow-Avoided",
                               .$time < time.mid & .$sel > 0.5 ~ "Fast-Preferred",
                               .$time < time.mid & .$sel < 0.5 ~ "Fast-Avoided"))



ggplot() +
  geom_raster(data = all.dat, aes(x, y, fill = fun_class), na.rm = T) +
  scale_fill_manual("", values = c("firebrick","forestgreen","yellow","steelblue1"),
                    na.translate = F) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.85, color = "black") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  facet_wrap(~ season, nrow = 2)

# ggsave("Giant Armadillo Connectivity_disc.png", width = 7, height = 5,
#        units = "in", dpi = 300)
