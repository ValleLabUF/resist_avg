

### Plot Velocity-based Resistance Vs Habitat Selection

library(ggnewscale)

#############
### North ###
#############

# Rescale Habitat Selection on -1 to 1 scale (from 0-1 scale)

# hab.selN<- (ssfSurfN_s.df$sel - 0.5)/0.5
hab.selN<- ssfSurfN_s.df$sel

dat.compN<- data.frame(sel = hab.selN, res = resist.mean.N.df$time,
                       lulc = factor(raster::values(lulcN)),
                       ndvi = raster::values(ndviN))
levels(dat.compN$lulc)<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")
time.midN<- (max(dat.compN$res, na.rm = T) + min(dat.compN$res, na.rm = T))/2

ggplot(dat.compN, aes(res, sel, color = lulc)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.midN)) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "North Pantanal") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

# ggsave("N Pantanal 3-banded Armadillo Quadrant_lulc.png", width = 7, height = 5,
#        units = "in", dpi = 300)


ggplot(dat.compN, aes(res, sel, color = ndvi)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.midN)) +
  scale_color_distiller(palette = "RdBu", limits = c(0,1)) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "North Pantanal") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

# ggsave("N Pantanal 3-banded Armadillo Quadrant_ndvi.png", width = 7, height = 5,
#        units = "in", dpi = 300)



### Reclassify pixels and map

## Using 2 continuous color scales
ggplot() +
  geom_raster(data = resist.mean.N.df, aes(x, y, fill = time), alpha = 0.6) +
  scale_fill_distiller("Time Spent\nper Cell (min)", palette = "Blues",
                       na.value = "transparent", direction = 1) +
  new_scale_fill() +
  geom_raster(data = ssfSurfN_s.df, aes(x, y, fill = sel), alpha = 0.6) +
  scale_fill_distiller("Selection", palette = "Reds",
                       na.value = "transparent", direction = 1, limits = c(0,1)) +
  geom_path(data = dat %>% filter(region == "N"), aes(x, y, group = id), alpha = 0.5,
            color = "chartreuse") +
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
        legend.text = element_text(size = 12))

# ggsave("N Pantanal 3-banded Armadillo Connectivity_cont.png", width = 7, height = 5,
#        units = "in", dpi = 300)




## Creation of 4 classes

all.datN<- cbind(resist.mean.N.df, sel = ssfSurfN_s.df$sel)

all.datN<- all.datN %>% 
  mutate(fun_class = case_when(.$time > time.midN & .$sel > 0.5 ~ "Slow Connectivity",
                               .$time > time.midN & .$sel < 0.5 ~ "Impediment",
                               .$time < time.midN & .$sel > 0.5 ~ "Fast Connectivity",
                               .$time < time.midN & .$sel < 0.5 ~ "Risky"))



ggplot() +
  geom_raster(data = all.datN, aes(x, y, fill = fun_class)) +
  scale_fill_brewer("", palette = "Dark2", na.value = "transparent") +
  geom_path(data = dat %>% filter(region == "N"), aes(x, y, group = id), alpha = 0.35,
            color = "chartreuse") +
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
        legend.text = element_text(size = 12))

# ggsave("N Pantanal 3-banded Armadillo Connectivity_disc.png", width = 7, height = 5,
#        units = "in", dpi = 300)






#############
### South ###
#############

# hab.selS<- (ssfSurfS_s.df$sel - 0.5)/0.5
hab.selS<- ssfSurfS_s.df$sel

dat.compS<- data.frame(sel = hab.selS, res = resist.mean.S.df$time,
                       lulc = factor(raster::values(lulcS)),
                       ndvi = raster::values(ndviS))
levels(dat.compS$lulc)<- c("Field", "Forest", "Water", "Pasture", "Road")
time.midS<- (max(dat.compS$res, na.rm = T) + min(dat.compS$res, na.rm = T))/2

ggplot(dat.compS, aes(res, sel, color = lulc)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.midS)) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "South Pantanal") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

# ggsave("S Pantanal 3-banded Armadillo Quadrant_lulc.png", width = 7, height = 5,
#        units = "in", dpi = 300)


ggplot(dat.compS, aes(res, sel, color = ndvi)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.midS)) +
  scale_color_distiller(palette = "RdBu", limits = c(0,1)) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "South Pantanal") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

# ggsave("S Pantanal 3-banded Armadillo Quadrant_ndvi.png", width = 7, height = 5,
#        units = "in", dpi = 300)



### Reclassify pixels and map

## Using 2 continuous color scales
ggplot() +
  geom_raster(data = resist.mean.S.df, aes(x, y, fill = time), alpha = 0.6) +
  scale_fill_distiller("Time Spent\nper Cell (min)", palette = "Blues",
                       na.value = "transparent", direction = 1) +
  new_scale_fill() +
  geom_raster(data = ssfSurfS_s.df, aes(x, y, fill = sel), alpha = 0.6) +
  scale_fill_distiller("Selection", palette = "Reds",
                       na.value = "transparent", direction = 1, limits = c(0,1)) +
  geom_path(data = dat %>% filter(region == "S"), aes(x, y, group = id), alpha = 0.5,
            color = "chartreuse") +
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
        legend.text = element_text(size = 12))

# ggsave("S Pantanal 3-banded Armadillo Connectivity_cont.png", width = 7, height = 5,
#        units = "in", dpi = 300)




## Creation of 4 classes

all.datS<- cbind(resist.mean.S.df, sel = ssfSurfS_s.df$sel)

all.datS<- all.datS %>% 
  mutate(fun_class = case_when(.$time > time.midS & .$sel > 0.5 ~ "Slow Connectivity",
                               .$time > time.midS & .$sel < 0.5 ~ "Impediment",
                               .$time < time.midS & .$sel > 0.5 ~ "Fast Connectivity",
                               .$time < time.midS & .$sel < 0.5 ~ "Risky"))



ggplot() +
  geom_raster(data = all.datS, aes(x, y, fill = fun_class)) +
  scale_fill_brewer("", palette = "Dark2", na.value = "transparent") +
  geom_path(data = dat %>% filter(region == "S"), aes(x, y, group = id), alpha = 0.35,
            color = "chartreuse") +
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
        legend.text = element_text(size = 12))

# ggsave("S Pantanal 3-banded Armadillo Connectivity_disc.png", width = 7, height = 5,
#        units = "in", dpi = 300)