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
yan_conifer <- rast('outputs/spatial/Yan/yans_files_conifer.tif')
F3 <- rast('outputs/spatial/livetrees/TPA24.tif')
F3_conifer <- rast('outputs/spatial/livetrees/TPA24_conifer.tif')
salo <- rast('outputs/spatial/livetrees/salo_cover.tif')
salo_conifer <- rast('outputs/spatial/livetrees/salo_cover_conifer.tif')
climate <- rast('outputs/spatial/weather/climate_means_2000_2020_MayOct.tif')

# conifer = 1, hardwood = 2, barren = 3, shrub = 4, wetland = 5, water = 6, herb = 7, urb = 8, subalpine = 9
# other(not really burnable) = 3,5,6,8, conifer = 1,9, hardwood = 2, shrub/herb = 4,7
veg <- rast('data/SEKI_veg/veg30m.tif') 
rclmat <- matrix(
  c(0, NA,
    1, 1,
    2, 2, 
    4, 3, 
    7, 3,
    3, 4,
    5, 4,
    6, 4,
    8, 4,
    9, 1), 
  ncol = 2, byrow = T
)
veg_rc <- classify(veg, rclmat)
# now: 1 = conifer, 2 = hardwood, 3 = shrub/herbaceous, 4 = not really burnable 
# reproject any rasters ---------------------------------------------------

# climate needs to be resampled to a 30m grid. default is bilinear.
climate_crs <- project(climate, yan)



# last minute small manipulations -----------------------------------------

# change NAs to 0 in TPA and cover
F3 <- ifel(is.na(F3), 0, F3)
F3_conifer <- ifel(is.na(F3_conifer), 0, F3_conifer)
salo <- ifel(is.na(salo), 0, salo)
salo_conifer <- ifel(is.na(salo_conifer), 0, salo_conifer)


# merge -------------------------------------------------------------------

# gonna keep the two fires seperate. should lead to smaller overall space. 

# KNP ==================================

raster_list <- c(
  list(
    cbi_knp %>% setNames('cbi'),
    crop(firehistory_cbi, ndvi_knp) %>% setNames('FH_cbi'), 
    crop(firehistory_year, ndvi_knp) %>% setNames('FH_year'),
    ndvi_knp %>% setNames('ndvi'),
    crop(F3, ndvi_knp) %>% setNames('TPA'),
    crop(F3_conifer, ndvi_knp) %>% setNames('TPA_conifer'),
    crop(salo, ndvi_knp) %>% setNames('cover'),
    crop(salo_conifer, ndvi_knp) %>% setNames('cover_conifer')
  ) ,
  topo_knp,
  crop(yan, ndvi_knp) %>% setNames( names(.)),
  crop(yan_conifer, ndvi_knp) %>% setNames(glue::glue("{names(.)}_con")),
  crop(climate_crs, ndvi_knp) %>% setNames( names(.)), 
  crop(veg_rc, ndvi_knp) %>% setNames('veg')
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
    crop(F3, ndvi_castle) %>% setNames('TPA'),
    crop(F3_conifer, ndvi_castle) %>% setNames('TPA_conifer'),
    crop(salo, ndvi_castle) %>% setNames('cover'),
    crop(salo_conifer, ndvi_castle) %>% setNames('cover_conifer')
  ) ,
  topo_castle,
  crop(yan, ndvi_castle) %>% setNames( names(.)),
  crop(yan_conifer, ndvi_castle) %>% setNames(glue::glue("{names(.)}_con")),
  crop(climate_crs, ndvi_castle) %>% setNames( names(.)), 
  crop(veg_rc, ndvi_castle) %>% setNames('veg')
)

# force rasters to common extent. use topo layer.
rasters_castle <- map(raster_list, extend, y = ndvi_castle) %>% rast 



# export ------------------------------------------------------------------

path <- 'outputs/spatial/compiled'
if(!dir.exists(path)) dir.create(path, recursive = T)
writeRaster(rasters_knp, filename = file.path(path, 'rasters_knp.tif'), overwrite = T)
writeRaster(rasters_castle, filename = file.path(path, 'rasters_castle.tif'), overwrite = T)

