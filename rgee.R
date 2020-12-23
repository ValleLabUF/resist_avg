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

ee_Initialize()



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






# downConfig = list(scale = 30, maxPixels = 1.0E13, driveFolder = 'image')
# img_lst = naip_2015$toList(10)
# 
# #Method - 02
# usda_stars_time <- list()
# for (index in seq_len(count)) {
#   image = ee$Image(img_lst$get(index-1))
#   name = image$get('system:index')$getInfo()
#   # print(name)
#   usda_stars <- ee_as_stars(
#     image = image,
#     region = polys,
#     scale = downConfig$scale,
#     geodesic = FALSE,
#     maxPixels = downConfig$maxPixels,
#     container = downConfig$driveFolder
#   )
#   names(usda_stars) <- name
#   usda_stars[usda_stars==0] = NA
#   usda_stars_time[[index]] <- usda_stars
# }
# usda_stars_mosaic <- do.call(st_mosaic, usda_stars_time)




#load files as raster brick
ndvi<- list.files(getwd(), pattern = "*.tif$")
ndvi.stack<- stack(ndvi)
ndvi.brick<- brick(ndvi.stack)

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

mosaic.rast<- list()
tmp1<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndvi.stack))
ind.m<- format(as.Date(tmp1, format = "%Y%m%d"), format = "%m") %>% 
  as.numeric()
ord.months<- c(5:12,1)


for (i in 1:length(unique(ind.m))) {
  ind<- which(ind.m %in% ord.months[i])  #using unique(ind.m) instead of ord.months results in wrong order
  ndvi.sub<- ndvi.stack[[ind]]
  
  # time1<- gsub(pattern = "LC08_22.07._", replacement = "", names(ndvi.stack))[i]
  ind.pr<- stringr::str_sub(names(ndvi.stack)[ind], 6, 11)

  if (length(unique(ind.pr)) > 1) {
    tile.list<- vector("list", length = length(unique(ind.pr)))
    
    for (j in 1:length(unique(ind.pr))) {
      ind.tiles<- grep(pattern = unique(ind.pr)[j], names(ndvi.stack)[ind])
      tile.list[[j]]<- mean(ndvi.sub[[ind.tiles]], na.rm=TRUE)
    }
    
    tile.list$fun <- mean
    tile.list$na.rm <- TRUE
    
    tile.mosaic<- do.call(raster::mosaic, tile.list)
    
  } else {
    tile.mosaic<- mean(ndvi.sub, na.rm=TRUE)
  }
  
  mosaic.rast[[i]]<- tile.mosaic
  }


coordinates(dat) <- ~x+y


mosaic.rast<- raster::brick(mosaic.rast)
names(mosaic.rast)<- month.abb[c(5:12,1)]
rasterVis::levelplot(mosaic.rast, at=breaks, col.regions=cols, main="NDVI") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("red", 0.5)))



### If wanting to aggregate into seasons:

# May is Fall
# Jun/Jul/Aug are Winter
# Sep/Oct/Nov are Spring
# Dec/Jan are Summer

season.ind<- c("Fall", rep("Winter",3), rep("Spring",3), rep("Summer",2))
ndvi.season<- stackApply(mosaic.rast, season.ind, fun = mean)
names(ndvi.season)<- c("Fall","Winter","Spring","Summer")
rasterVis::levelplot(ndvi.season, at=breaks, col.regions=cols, main="NDVI") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("red", 0.5)))




### Try plotting in ggplot2 so points can be plotted by month and season

mosaic.rast.df<- as.data.frame(mosaic.rast, xy=TRUE)
mosaic.rast.df<- pivot_longer(mosaic.rast.df, cols = -c(x,y), names_to = "month",
                              values_to = "ndvi")
mosaic.rast.df$month<- factor(mosaic.rast.df$month, levels = month.abb[c(5:12,1)])
mosaic.rast.df<- filter(mosaic.rast.df, ndvi > 0)  #only keep values from 0-1

dat<- as.data.frame(dat)
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))


ggplot() +
  geom_raster(data = mosaic.rast.df, aes(x, y, fill = ndvi)) +
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




####################################################
### Write and save RasterBricks as geoTIFF files ###
####################################################

setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

## Monthly

#if saving individual layers
# writeRaster(mosaic.rast, filename=names(mosaic.rast), bylayer=TRUE, format="GTiff")

#if saving RasterBrick as raster format
writeRaster(mosaic.rast, filename = 'GiantArm_ndvi_monthly.grd', format="raster",
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

#if saving RasterBrick as geoTIFF
# writeRaster(ndvi.season, filename='GiantArm_ndvi_season.tif', format="GTiff", overwrite=TRUE,
#             options=c("INTERLEAVE=BAND","COMPRESS=LZW"))