# creating a MMI layer reporting events from Aug 2015 - July 2016


library(sf)
library(tidyverse)
library(raster)
library(lubridate)
theme_set(theme_classic())

evty15 <- raster('../../../Documents/ArcGIS/Spatial_data/edart/evty2015_northbox.tif')
evty16 <- raster('../../../Documents/ArcGIS/Spatial_data/edart/evty2016_northbox.tif')
evty17 <- raster('../../../Documents/ArcGIS/Spatial_data/edart/evty2017_northbox.tif')
mmi15 <- raster('../../../Documents/ArcGIS/Spatial_data/edart/mmi2015_northbox.tif')
mmi16 <- raster('../../../Documents/ArcGIS/Spatial_data/edart/mmi2016_northbox.tif')
mmi17 <- raster('../../../Documents/ArcGIS/Spatial_data/edart/mmi2017_northbox.tif')

# Flight dates: Aug 2015 - July 2016. TAOs are mortality in that interval. 
# Match MMI mortality events for those dates. 

sum(values(evty15 == -1)) #858 pixs need to be removed
sum(values(evty16 == -1)) #858 pixs need to be removed
sum(values(evty17 == -1)) #858 pixs need to be removed

# MMI 2015 ----------------------------------------------------------------

# filter MMI 2015 to include events Aug 1 - Dec 31
unique(values(evty15)) %>% sort %>% date_decimal()
#[1]   -1.000    0.000 2015.360 2015.403 2015.447 2015.491 2015.535 2015.579 2015.622
#[10] 2015.666 2015.710 2015.798 2015.885

# convert to decimal dates
start15 <- ymd('2015-08-01') %>% decimal_date()

#mask out pixels before start date and -1. Keep non-events. Apply mask to MMI. 
mask15 <- evty15 == 0 | evty15 > start15 #things to keep
mmi15_m <- mask(mmi15, mask15, maskvalue = 0)
plot(mmi15_m)

#include everything from 2015 except pixels to be ignored
mask15v2 <- evty15 >= 0 #things to keep
mmi15_mv2 <- mask(mmi15, mask15v2, maskvalue = 0)
plot(mmi15_mv2)

sum(values(mmi15_m == 0), na.rm = T) #871538 zero values
sum(values(mmi15_mv2 == 0), na.rm = T) 


# MMI 2016 ----------------------------------------------------------------

# create two versions:
# 1) filter MMI 2016 to include events Jan 1 - July 31
# 2) No filter--last events recorded so filtering it will possibly underestimate mortality
# unique(values(evty16))
# [1]    0.000 2016.367 2016.499 2016.849 2016.630 2016.455 2016.586
# [8] 2016.674 2016.411 2016.718 2016.542 2016.805   -1.000

# convert to decimal dates
end16 <- ymd('2016-07-31') %>% decimal_date()

#mask out pixels not in range. Apply mask to MMI. 
mask16_v1 <- evty16 >= 0L & evty16 <= end16 #things to keep
mmi16_v1 <- mask(mmi16, mask16_v1, maskvalue = 0)
mask16_v2 <- evty16 >= 0L #things to keep--all pixs except for ignored
mmi16_v2 <- mask(mmi16, mask16_v2, maskvalue = 0)

sum(values(mmi16_v2 == 0), na.rm = T) #836737 zero values
plot(mmi16_v2)

# MMI 2017 ----------------------------------------------------------------

# also get the first 2 images of 2017 in case mortality was missed in 2016 from smoke
dates17 <- unique(values(evty17)) %>% sort
#[1]   -1.000    0.000 2017.462 2017.506 2017.550 2017.638 2017.725 2017.769 2017.813
#[10] 2017.857 2017.900
mask17 <- evty17 == 0L | evty17 %in% dates17[3:4] #things to keep
mmi17_v1 <- mask(mmi17, mask17, maskvalue = 0)
plot(mmi17_v1)


# MERGE 2015 and 2016 MMI ---------------------------------------------------

# merge by summing pixels? 
MMI_v1 <- sum(mmi15_m, mmi16_v1, na.rm = F) #includes late 2015, early 2016
MMI_v2 <- sum(mmi15_m, mmi16_v2, na.rm = F) #includes late 2015, all of 2016
MMI_v22 <- sum(mmi15_mv2, mmi16_v2, na.rm = F) #includes all of 2015 and 2016
MMI_v221 <- sum(mmi15_mv2, mmi16_v2, na.rm = F) #includes all of 2015, 2016, and first 2 images of 2017

# truncate values at 100
MMI_v22[MMI_v22 > 100] <- 100
MMI_v2[MMI_v2 > 100] <- 100
MMI_v221[MMI_v221 > 100] <- 100

plot(MMI_v2)
plot(MMI_v22)
plot(MMI_v221) #hardly any different from v22

names(MMI_v1) <- 'MMI_v1'
names(MMI_v2) <- 'MMI_v2'
names(MMI_v22) <- 'MMI_v22'
names(MMI_v221) <- 'MMI_v221'
sum(values(MMI_v2 == 0), na.rm = T) #811363 zero values. includes pixels with non-events. 

hist(MMI_v22$MMI_v22[MMI_v22>0])
hist(MMI_v22) #that includes values outside the concave hull and non-events inside it.

# export mmi raster
writeRaster(MMI_v1, filename = "outputs/rasters/MMI_v1.tif", overwrite = T)
writeRaster(MMI_v2, filename = "outputs/rasters/MMI_v2.tif", overwrite = T)
writeRaster(MMI_v22, filename = "outputs/rasters/MMI_v22.tif", overwrite = T)
writeRaster(MMI_v221, filename = "outputs/rasters/MMI_v221.tif", overwrite = T)

# CREATE COMMON EVTY RASTER TOO -------------------------------------------

# I want to know when the pixels were detected and if a pixel has 2 events, use latest one
evty15m <- mask(evty15, mask15v2, maskvalue = 0)
evty16m <- mask(evty16, mask16_v2, maskvalue = 0)
evty17m <- mask(evty17, mask17, maskvalue = 0)
evtycommon <- stack(evty15m, evty16m, evty17m) %>% max
writeRaster(evtycommon, filename = "outputs/rasters/evtycommon_v221.tif", overwrite = T)
