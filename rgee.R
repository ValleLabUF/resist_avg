## Landsat8 collected on 16-day interval


library(tidyverse)
library(sf)
library(raster)
library(lubridate)
library(sp)
library(furrr)
library(future)
library(rgee)
library(raster)
library(rasterVis)
library(googledrive)

ee_Initialize(email = "joshcullen10@gmail.com")



setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")
dat<- read.csv("Binned Armadillo Acceleration Data.csv", as.is = T)
dat$date<- as_datetime(dat$date)

# Filter out observations where coords are NA
dat<- dat %>% 
  filter(!is.na(x))


## Date range for giant armadillos
range(dat$date)  #2019-05-19 to 2020-01-24


setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")
tmp<- raster('cheia_UTM1.tif')
tmp1<- raster('cheiann_UTM1.tif')



# Create dummy rasters for extent of each region (with buffer)

#North
rast<- tmp
values(rast)<- 0



#################################################################
### Download data from Google Earth Engine for North Pantanal ###
#################################################################

# setwd("~/Documents/Snail Kite Project/Data/armadillos/NDVI/North Pantanal")
# 
# ndvi.N<- dir(getwd(), "*.tif$")
# for (i in ndvi.N) assign(i, raster(i))



bounds<- st_as_sf(data.frame(rasterToPoints(rast)), coords = c("x","y"), 
                    crs = "+init=epsg:32721") %>%  #convert to sf object
  st_bbox() %>%  #extract bounding box
  st_as_sfc() %>%  #convert into polygon
  sf_as_ee()  #convert into GEE format 



# This function gets NDVI from Landsat 8 imagery
addNDVI<- function(image) {
  return(image$addBands(image$normalizedDifference(c("B5", "B4"))$
                          rename('NDVI')))
}

# This function gets NDWI from Landsat 8 imagery
addNDWI<- function(image) {
  return(image$addBands(image$normalizedDifference(c("B5", "B6"))$
                          rename('NDWI')))
}

# This function calculates AWEIsh from Landsat 8 imagery
addAWEIsh<- function(image) {
  
  return(image$addBands(image$expression(
    expression = '(BLUE + (2.5 * GREEN)) - (1.5 * (NIR + SWIR1)) - (0.25 * SWIR2)',
    opt_map =  list(
      'SWIR2' = image$select('B7'),
      'SWIR1' = image$select('B6'),
      'NIR' = image$select('B5'),
      'GREEN' = image$select('B3'),
      'BLUE' = image$select('B2')
    )
  )$rename('AWEI')))
  
}

# This functions masks clouds and cloud shadows (bits 3 and 5)
maskL8sr<- function(image) {
  cloudShadowBitMask<- bitwShiftL(1, 3)
  cloudsBitMask<- bitwShiftL(1, 5)
  qa<- image$select('pixel_qa')  # Get the pixel QA band
  # Both flags should be set to zero, indicating clear conditions
  mask<- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$
    And(qa$bitwiseAnd(cloudsBitMask)$eq(0))
  return(image$updateMask(mask))
}



# Retrieve Landsat 8-based NDVI
l8<- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
  filterBounds(bounds)$
  filterDate('2019-05-19', '2020-01-24')$
  map(function(image){image$clip(bounds)})$
  sort("system:time_start", TRUE)$  # Sort the collection in chronological order
  map(function(x) x$reproject("EPSG:32721"))

print(l8$size()$getInfo())  #check number of images




# Map the NDVI and cloud mask functions; select only the NDVI band
l8_ndvi<- l8$map(addNDVI)$
  map(maskL8sr)$
  select('NDVI')

ee_print(l8_ndvi)


# Plot the map of the median NDVI for the region
Map$setCenter(-55.76664, -19.19482, 11)
Map$addLayer(l8_ndvi$median(),
             visParams = list(
               min = 0.0,
               max = 1.0,
               bands = "NDVI",
               palette = c(
                 'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
                 '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
                 '012E01', '011D01', '011301'
               )
             ))




# Map the AWEI and cloud mask functions; select only the AWEI band
# l8_awei<- l8$map(addAWEIsh)$
#   map(maskL8sr)$
#   select('AWEI')
# 
# ee_print(l8_awei)
# mean1<- l8_awei$first()$reduceRegions(collection =  bounds,
#                               reducer = ee$Reducer$mean(),
#                               scale = 30)
# print(mean1$getInfo())
# 
# ndwiViz <- list(min = 0, max = 1, palette = c("00FFFF", "0000FF"))
# 
# # Plot the map of the median AWEI for the region
# Map$setCenter(-55.76664, -19.19482, 11)
# Map$addLayer(l8_awei$mean(),
#              visParams = ndwiViz)


# Map the NDWI and cloud mask functions; select only the NDWI band
l8_ndwi<- l8$map(addNDWI)$
  map(maskL8sr)$
  select('NDWI')

ee_print(l8_ndwi)
print(l8_ndwi$first()$getInfo())

ndwiViz <- list(min = 0, max = 1, palette = c("00FFFF", "0000FF"))

# Plot the map of the median NDWI for the region
Map$setCenter(-55.76664, -19.19482, 11)
Map$addLayer(l8_ndwi$median(),
             visParams = ndwiViz)







###########################
### Export Data Locally ###
###########################

# Download all images locally
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg/NDVI")

ndvi<- ee_imagecollection_to_local(
  ic = l8_ndvi,
  scale = 30,
  region = bounds,
  via = 'drive'
)


ndvi.stack<- raster::stack(ndvi)

breaks<- seq(0, 1, by=0.01)
cols<- colorRampPalette(c('#FFFFFF', '#CE7E45', '#DF923D', '#F1B555', '#FCD163', '#99B718',
                          '#74A901', '#66A000', '#529400', '#3E8601', '#207401', '#056201',
                          '#004C00', '#023B01', '#012E01', '#011D01', '#011301'
))(length(breaks)-1)

##plot
rasterVis::levelplot(ndvi.stack, at=breaks, col.regions=cols, main="NDVI")









# # Download all images locally
# setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg/AWEI")
# 
# awei<- ee_imagecollection_to_local(
#   ic = l8_awei,
#   scale = 30,
#   region = bounds,
#   via = 'drive'
# )
# 
# 
# awei.stack<- raster::stack(awei)
# 
# breaks<- seq(0, 1, by=0.01)
# cols<- colorRampPalette(c("#00FFFF", "#0000FF"))(length(breaks)-1)
# 
# ##plot
# rasterVis::levelplot(awei.stack[[1:4]], at=breaks, col.regions=cols, main="AWEI")






# Download all images locally
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg/NDWI")

ndwi<- ee_imagecollection_to_local(
  ic = l8_ndwi,
  scale = 30,
  region = bounds,
  via = 'drive'
)


ndwi.stack<- raster::stack(ndwi)
ndwi.stack<- raster::flip(ndwi.stack, direction = 'y')

breaks<- seq(0, 1, by=0.01)
cols<- colorRampPalette(c("#00FFFF", "#0000FF"))(length(breaks)-1)

##plot
rasterVis::levelplot(ndwi.stack[[1:4]], at=breaks, col.regions=cols, main="NDWI")






###############################################
### Mosaic and Average Data by Month/Season ###
###############################################


#### NDVI ####

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg/NDVI")

#load files as raster brick
ndvi<- list.files(getwd(), pattern = "*.tif$")
ndvi.stack<- stack(ndvi)
ndvi.stack2<- flip(ndvi.stack, direction = 'y')
names(ndvi.stack2)<- names(ndvi.stack)
ndvi.brick<- ndvi.stack2

#change values x == 0 to NA (these are masked pixels)
raster::values(ndvi.brick)[raster::values(ndvi.brick) == 0] <- NA

#rename images by date
names(ndvi.brick)<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndvi.brick))


#visualize original rasters
breaks<- seq(0, 1, by=0.01)
cols<- colorRampPalette(c('#FFFFFF', '#CE7E45', '#DF923D', '#F1B555', '#FCD163', '#99B718',
                          '#74A901', '#66A000', '#529400', '#3E8601', '#207401', '#056201',
                          '#004C00', '#023B01', '#012E01', '#011D01', '#011301'
))(length(breaks)-1)
rasterVis::levelplot(ndvi.brick[[1:6]], at=breaks, col.regions=cols, main="NDVI")






#### Mosaic RasterBrick over space and time ####

mosaic.ndvi<- list()
tmp1<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndvi.stack2))
ind.m<- format(as.Date(tmp1, format = "%Y%m%d"), format = "%m") %>% 
  as.numeric()
ord.months<- c(5:12,1)


for (i in 1:length(unique(ind.m))) {
  ind<- which(ind.m %in% ord.months[i])  #using unique(ind.m) instead of ord.months results in wrong order
  ndvi.sub<- ndvi.stack2[[ind]]
  
  # time1<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndvi.stack))[i]
  ind.pr<- stringr::str_sub(names(ndvi.stack2)[ind], 6, 11)

  if (length(unique(ind.pr)) > 1) {
    tile.list<- vector("list", length = length(unique(ind.pr)))
    
    for (j in 1:length(unique(ind.pr))) {
      ind.tiles<- grep(pattern = unique(ind.pr)[j], names(ndvi.stack2)[ind])
      tile.list[[j]]<- mean(ndvi.sub[[ind.tiles]], na.rm=TRUE)
    }
    
    tile.list$fun <- mean
    tile.list$na.rm <- TRUE
    
    tile.mosaic<- do.call(raster::mosaic, tile.list)
    
  } else {
    tile.mosaic<- mean(ndvi.sub, na.rm=TRUE)
  }
  
  mosaic.ndvi[[i]]<- tile.mosaic
  }


coordinates(dat) <- ~x+y


mosaic.ndvi<- raster::brick(mosaic.ndvi)
names(mosaic.ndvi)<- month.abb[c(5:12,1)]
rasterVis::levelplot(mosaic.ndvi, at=breaks, col.regions=cols, main="NDVI") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("red", 0.5)))



### If wanting to aggregate into seasons:

# May is Fall
# Jun/Jul/Aug are Winter
# Sep/Oct/Nov are Spring
# Dec/Jan are Summer

season.ind<- c("Fall", rep("Winter",3), rep("Spring",3), rep("Summer",2))
ndvi.season<- stackApply(mosaic.ndvi, season.ind, fun = mean)
names(ndvi.season)<- c("Fall","Winter","Spring","Summer")
rasterVis::levelplot(ndvi.season, at=breaks, col.regions=cols, main="NDVI") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("red", 0.5)))




### Try plotting in ggplot2 so points can be plotted by month and season

mosaic.ndvi.df<- as.data.frame(mosaic.ndvi, xy=TRUE)
mosaic.ndvi.df<- pivot_longer(mosaic.ndvi.df, cols = -c(x,y), names_to = "month",
                              values_to = "ndvi")
mosaic.ndvi.df$month<- factor(mosaic.ndvi.df$month, levels = month.abb[c(5:12,1)])
mosaic.ndvi.df<- filter(mosaic.ndvi.df, ndvi > 0)  #only keep values from 0-1

dat<- as.data.frame(dat)
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))


ggplot() +
  geom_raster(data = mosaic.ndvi.df, aes(x, y, fill = ndvi)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, values = breaks, na.value = "transparent") +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ month) +
  theme(strip.text = element_text(size = 12, face = "bold"))




ndvi.season.df<- as.data.frame(ndvi.season, xy=TRUE)
ndvi.season.df<- pivot_longer(ndvi.season.df, cols = -c(x,y), names_to = "season",
                              values_to = "ndvi")
ndvi.season.df$season<- factor(ndvi.season.df$season,
                               levels = c("Fall","Winter","Spring","Summer"))
ndvi.season.df<- filter(ndvi.season.df, ndvi > 0)  #only keep values from 0-1



ggplot() +
  geom_raster(data = ndvi.season.df, aes(x, y, fill = ndvi)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, values = breaks, na.value = "transparent") +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ season) +
  theme(strip.text = element_text(size = 12, face = "bold"))







#### NDWI ####

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg/NDWI")

#load files as raster brick
ndwi<- list.files(getwd(), pattern = "*.tif$")
ndwi.stack<- stack(ndwi)
ndwi.stack2<- flip(ndwi.stack, direction = 'y')
names(ndwi.stack2)<- names(ndwi.stack)
ndwi.brick<- ndwi.stack2

#change values x == 0 to NA (these are masked pixels)
raster::values(ndwi.brick)[raster::values(ndwi.brick) == 0] <- NA

#rename images by date
names(ndwi.brick)<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndwi.brick))


#visualize original rasters (modified Zissou1 palette from {wesanderson})
breaks<- seq(-1, 1, by=0.01)
cols<- colorRampPalette(rev(c("#3B9AB2","#EBCC2A","#F21A00")))(length(
  breaks)-1)
rasterVis::levelplot(ndwi.brick[[1:6]], at=breaks, col.regions=cols, main="NDWI")




#### Mosaic RasterBrick over space and time ####

mosaic.ndwi<- list()
tmp1<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndwi.stack2))
ind.m<- format(as.Date(tmp1, format = "%Y%m%d"), format = "%m") %>% 
  as.numeric()
ord.months<- c(5:12,1)


for (i in 1:length(unique(ind.m))) {
  ind<- which(ind.m %in% ord.months[i])  #using unique(ind.m) instead of ord.months results in wrong order
  ndwi.sub<- ndwi.stack2[[ind]]
  
  # time1<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndvi.stack))[i]
  ind.pr<- stringr::str_sub(names(ndwi.stack2)[ind], 6, 11)
  
  if (length(unique(ind.pr)) > 1) {
    tile.list<- vector("list", length = length(unique(ind.pr)))
    
    for (j in 1:length(unique(ind.pr))) {
      ind.tiles<- grep(pattern = unique(ind.pr)[j], names(ndwi.stack2)[ind])
      tile.list[[j]]<- mean(ndwi.sub[[ind.tiles]], na.rm=TRUE)
    }
    
    tile.list$fun <- mean
    tile.list$na.rm <- TRUE
    
    tile.mosaic<- do.call(raster::mosaic, tile.list)
    
  } else {
    tile.mosaic<- mean(ndwi.sub, na.rm=TRUE)
  }
  
  mosaic.ndwi[[i]]<- tile.mosaic
}


coordinates(dat) <- ~x+y


mosaic.ndwi<- raster::brick(mosaic.ndwi)
names(mosaic.ndwi)<- month.abb[c(5:12,1)]
rasterVis::levelplot(mosaic.ndwi, at=breaks, col.regions=cols, main="NDWI") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("red", 0.5)))



### If wanting to aggregate into seasons:

# May is Fall
# Jun/Jul/Aug are Winter
# Sep/Oct/Nov are Spring
# Dec/Jan are Summer

season.ind<- c("Fall", rep("Winter",3), rep("Spring",3), rep("Summer",2))
ndwi.season<- stackApply(mosaic.ndwi, season.ind, fun = mean)
names(ndwi.season)<- c("Fall","Winter","Spring","Summer")
rasterVis::levelplot(ndwi.season, at=breaks, col.regions=cols, main="NDWI") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("red", 0.5)))




### Try plotting in ggplot2 so points can be plotted by month and season

mosaic.ndwi.df<- as.data.frame(mosaic.ndwi, xy=TRUE)
mosaic.ndwi.df<- pivot_longer(mosaic.ndwi.df, cols = -c(x,y), names_to = "month",
                              values_to = "ndwi")
mosaic.ndwi.df$month<- factor(mosaic.ndwi.df$month, levels = month.abb[c(5:12,1)])
# mosaic.ndwi.df<- filter(mosaic.ndwi.df, ndwi > 0)  #only keep values from 0-1

dat<- as.data.frame(dat)
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))




ggplot() +
  geom_raster(data = mosaic.ndwi.df, aes(x, y, fill = ndwi)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(-1,1)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ month) +
  theme(strip.text = element_text(size = 12, face = "bold"))




ndwi.season.df<- as.data.frame(ndwi.season, xy=TRUE)
ndwi.season.df<- pivot_longer(ndwi.season.df, cols = -c(x,y), names_to = "season",
                              values_to = "ndwi")
ndwi.season.df$season<- factor(ndwi.season.df$season,
                               levels = c("Fall","Winter","Spring","Summer"))
# ndwi.season.df<- filter(ndwi.season.df, ndwi > 0)  #only keep values from 0-1



ggplot() +
  geom_raster(data = ndwi.season.df, aes(x, y, fill = ndwi)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(-1,1)) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~ season) +
  labs(x="Easting", y="Northing", title = "Pantanal NDWI") +
  coord_equal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# ggsave("Giant Armadillo NDWI_seasons.png", width = 6, height = 5, units = "in", dpi = 300)



####################################################
### Write and save RasterBricks as geoTIFF files ###
####################################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

## Monthly

#if saving individual layers
# writeRaster(mosaic.rast, filename=names(mosaic.rast), bylayer=TRUE, format="GTiff")

#if saving RasterBrick as raster format
writeRaster(mosaic.ndvi, filename = 'GiantArm_ndvi_monthly.grd', format="raster",
            overwrite=TRUE)
writeRaster(mosaic.ndwi, filename = 'GiantArm_ndwi_monthly.grd', format="raster",
            overwrite=TRUE)

#if saving RasterBrick as geoTIFF
# writeRaster(mosaic.rast, filename='GiantArm_ndvi_monthly.tif', format="GTiff", overwrite=TRUE,
#             options=c("INTERLEAVE=BAND","COMPRESS=LZW"))



## Season

#if saving individual layers
# writeRaster(ndvi.season, filename=names(ndvi.season), bylayer=TRUE, format="GTiff")

#if saving RasterBrick as raster format
writeRaster(ndvi.season, filename = 'GiantArm_ndvi_season.grd', format="raster",
            overwrite=TRUE)
writeRaster(ndwi.season, filename = 'GiantArm_ndwi_season.grd', format="raster",
            overwrite=TRUE)

#if saving RasterBrick as geoTIFF
# writeRaster(ndvi.season, filename='GiantArm_ndvi_season.tif', format="GTiff", overwrite=TRUE,
#             options=c("INTERLEAVE=BAND","COMPRESS=LZW"))