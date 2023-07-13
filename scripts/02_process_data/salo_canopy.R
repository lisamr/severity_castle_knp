# look at salo's canopy cover 

rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)

salo <- rast('data/Salo/TulareCounty-Vegetation-CanopyCover-2020-Summer-00010m.tif')



# reproject the data ------------------------------------------------------


# needs to be consistent with the rest of the rasters. use yans data.
yans_files <- rast('outputs/spatial/Yan/yans_files_conifer.tif')

salo_reproj <- project(salo, yans_files$dead_area)
names(salo_reproj) <- "live_cover"



# mask out non-conifers ---------------------------------------------------

# yan's models picked up a lot of senescing hardwoods. to be consistent, 
# mask out this too
veg <- rast('data/SEKI_veg/veg30m.tif') # confer = 1
conifer <- ifel(veg == 1, 1, 0)
conifer <- crop(conifer, salo_reproj)

# mask: if conifer==0, mask and update yan's files as 0. 
salo_conifer <- terra::mask(salo_reproj, conifer, 
                            maskvalue = 0, updatevalue = 0)



# export ------------------------------------------------------------------


if(!dir.exists('outputs/spatial/livetrees')) dir.create('outputs/spatial/livetrees')
writeRaster(salo_reproj, 'outputs/spatial/livetrees/salo_cover.tif')
writeRaster(salo_conifer, 'outputs/spatial/livetrees/salo_cover_conifer.tif')






# viz ---------------------------------------------------------------------
# 
par(mfrow = c(1,3))
plot(salo_conifer, main = 'salo')
plot(rast('outputs/spatial/livetrees/TPA24.tif'), main = 'F3')
plot(yans_files$dead_area, main = 'dead, yan')
