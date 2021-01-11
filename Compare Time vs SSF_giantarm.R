### Must first have run Step 3_Analyze Results_giantarm.R and giantarm_ssf.R already


library(ggnewscale)



### Plot Velocity-based Resistance vs Habitat Selection

dat.comp<- data.frame(sel = ssf.pop2$mu, res = resist.pop2$mu,
                      lulc = factor(raster::values(lulc)),
                      ndwi = raster::values(mean(ndwi, na.rm = T)))
levels(dat.comp$lulc)<- c("Forest", "Closed_Savanna", "Open_Savanna", "Floodable")

time.mid<- (max(dat.comp$res, na.rm = T) + min(dat.comp$res, na.rm = T))/2

ggplot(dat.comp, aes(res, sel, color = lulc)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.mid)) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "Giant Armadillo") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

# ggsave("Giant Armadillo Quadrant_lulc.png", width = 7, height = 5,
#        units = "in", dpi = 300)



ggplot(dat.comp, aes(res, sel, color = ndwi)) +
  geom_hline(aes(yintercept = 0.5)) +
  geom_vline(aes(xintercept = time.mid)) +
  scale_color_distiller(palette = "RdBu", limits = c(-1,1)) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "Giant Armadillo") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

# ggsave("Giant Armadillo Quadrant_ndwi.png", width = 7, height = 5,
#        units = "in", dpi = 300)





### Reclassify pixels and map


## Using 2 continuous color scales
ggplot() +
  geom_raster(data = resist.pop2, aes(x, y, fill = mu), alpha = 0.6) +
  scale_fill_distiller("Time Spent\nper Cell (min)", palette = "Blues",
                       na.value = "transparent", direction = 1) +
  new_scale_fill() +
  geom_raster(data = ssf.pop2, aes(x, y, fill = mu), alpha = 0.6) +
  scale_fill_distiller("Selection", palette = "Reds",
                       na.value = "transparent", direction = -1, limits = c(0,1)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.5, color = "chartreuse") +
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

# ggsave("Giant Armadillo Connectivity_cont.png", width = 7, height = 5,
#        units = "in", dpi = 300)




## Creation of 4 classes

all.dat<- cbind(resist.pop2[,-4], ssf.pop2$mu)
names(all.dat)[3:4]<- c("time", "sel")

all.dat<- all.dat %>% 
  mutate(fun_class = case_when(.$time > time.mid & .$sel > 0.5 ~ "Slow Connectivity",
                                .$time > time.mid & .$sel < 0.5 ~ "Impediment",
                                .$time < time.mid & .$sel > 0.5 ~ "Fast Connectivity",
                                .$time < time.mid & .$sel < 0.5 ~ "Risky"))



ggplot() +
  geom_raster(data = all.dat, aes(x, y, fill = fun_class)) +
  scale_fill_brewer("", palette = "Dark2", na.value = "transparent") +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.35, color = "chartreuse") +
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

# ggsave("Giant Armadillo Connectivity_disc.png", width = 7, height = 5,
#        units = "in", dpi = 300)
