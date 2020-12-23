

### Plot Velocity-based Resistance Vs Habitat Selection

# Rescale Habitat Selection on -1 to 1 scale (from 0-1 scale)

hab.selN<- (ssfSurfN_s.df$sel - 0.5)/0.5

dat.compN<- data.frame(sel = hab.selN, res = resist.mean.N.df$time, lulc = factor(lulcN.df$lulc))
levels(dat.compN$lulc)<- c("Pasture", "HQ", "Fence", "Water", "Cane", "Forest")


ggplot(dat.compN, aes(res, sel, color = lulc)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point() +
  theme_bw() +
  ylim(-1,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "North Pantanal") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))





hab.selS<- (ssfSurfS_s.df$sel - 0.5)/0.5

dat.compS<- data.frame(sel = hab.selS, res = resist.mean.S.df$time, lulc = factor(lulcS.df$lulc))
levels(dat.compS$lulc)<- c("Field", "Forest", "Water", "Pasture", "Road")


ggplot(dat.compS, aes(res, sel, color = lulc)) +
  geom_hline(aes(yintercept = 0)) +
  geom_point() +
  theme_bw() +
  ylim(-1,1) +
  labs(x= "Time (min)", y="Habitat Preference", title = "South Pantanal") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 10))

