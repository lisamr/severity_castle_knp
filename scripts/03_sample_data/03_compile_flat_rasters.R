# compile geospatial data into a raster stack
library(terra)
library(tidyverse)
library(sf)

# load raster data inputs -------------------------------------------------

cbi_castle <- rast('outputs/spatial/GEE_CBI/issues/Castle_2020_CBI_bc_corrected.tif')
cbi_knp <- rast('outputs/spatial/GEE_CBI/KNP_2021_CBI_bc_corrected.tif')
firehistory_cbi <- rast('outputs/spatial/GEE_CBI/fireHistory_1985_2020_CBI_bc_corrected.tif')
firehistory_year <- rast('outputs/spatial/GEE_CBI/fireHistory_1985_2020_Fire_Year.tif')
topo_castle <- rast('outputs/spatial/topography/castle_terrain.tif')
topo_knp <- rast('outputs/spatial/topography/knp_terrain.tif')
ndvi_castle <- rast('data/indices/CASTLE_ndvi.tif')
ndvi_knp <- rast('data/indices/KNP_Complex_ndvi.tif')
yan <- rast('outputs/spatial/Yan/yans_files.tif')
F3 <- rast('outputs/spatial/livetrees/TPA24.tif')
climate <- rast('outputs/spatial/weather/climate_means_2000_2020_MayOct.tif')



# reproject any rasters ---------------------------------------------------

# climate needs to be resampled to a 30m grid. default is bilinear.
climate_crs <- project(climate, yan)


# merge -------------------------------------------------------------------

# gonna keep the two fires seperate. should lead to smaller overall space. 

# KNP ==================================

raster_list <- c(
  list(
    cbi_knp %>% setNames('cbi'),
    crop(firehistory_cbi, ndvi_knp) %>% setNames('FH_cbi'), 
    crop(firehistory_year, ndvi_knp) %>% setNames('FH_year'),
    ndvi_knp %>% setNames('ndvi'),
    crop(F3, ndvi_knp) %>% setNames('TPA')
  ) ,
  topo_knp,
  crop(yan, ndvi_knp) %>% setNames( names(.)),
  crop(climate_crs, ndvi_knp) %>% setNames( names(.))
  )


# force rasters to common extent. use ndvi layer.
rasters_knp <- map(raster_list, extend, y = ndvi_knp) %>% rast 


# Castle ==============================
raster_list <- c(
  list(
    cbi_castle %>% setNames('cbi'),
    crop(firehistory_cbi, ndvi_castle) %>% setNames('FH_cbi'), 
    crop(firehistory_year, ndvi_castle) %>% setNames('FH_year'),
    ndvi_castle %>% setNames('ndvi'),
    crop(F3, ndvi_castle) %>% setNames('TPA')
  ) ,
  topo_castle,
  crop(yan, ndvi_castle) %>% setNames( names(.)),
  crop(climate_crs, ndvi_castle) %>% setNames( names(.))
)

map(raster_list, ext)

# force rasters to common extent. use topo layer.
rasters_castle <- map(raster_list, extend, y = ndvi_castle) %>% rast 



# export ------------------------------------------------------------------

path <- 'outputs/spatial/compiled'
if(!dir.exists(path)) dir.create(path, recursive = T)
writeRaster(rasters_knp, filename = file.path(path, 'rasters_knp.tif'), overwrite = T)
writeRaster(rasters_castle, filename = file.path(path, 'rasters_castle.tif'), overwrite = T)

