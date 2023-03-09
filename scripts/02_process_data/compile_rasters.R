# compile geospatial data into a raster stack
library(terra)
library(tidyverse)
library(sf)

# load raster data inputs -------------------------------------------------

cbi_castle <- rast('outputs/spatial/GEE_CBI/Castle_2020_CBI_bc.tif')
cbi_knp <- rast('outputs/spatial/GEE_CBI/KNP_2021_CBI_bc.tif')
mtbs <- rast('data/fires/mtbs_SQF_KNP.grd')
firehistory_cbi <- rast('outputs/spatial/GEE_CBI/fireHistory_1985_2020_CBI_bc.tif')
firehistory_year <- rast('outputs/spatial/GEE_CBI/fireHistory_1985_2020_Fire_Year.tif')
firDir_castle <- rast('outputs/spatial/fire_direction/fire_dir_castle.tif')
firDir_knp <- rast('outputs/spatial/fire_direction/fire_dir_knp.tif')
topo_castle <- rast('outputs/spatial/topography/castle_terrain.tif')
topo_knp <- rast('outputs/spatial/topography/knp_terrain.tif')
ndvi_castle <- rast('data/indices/CASTLE_ndvi.tif')
ndvi_knp <- rast('data/indices/KNP_Complex_ndvi.tif')


# rasterize IR ------------------------------------------------------------

# load IR
IR_castle <- st_read('outputs/spatial/IR/castle_byDay_u.geojson')
IR_knp <- st_read('outputs/spatial/IR/knp_byDay_u.geojson')


f_rasterize_IR <- function(IR_data, topo_data){
  IR_data_r <- rasterize(rev(vect(IR_data)), topo_data, 'jday') 
  return(IR_data_r)
}

IR_rasters <- map2(list(IR_castle, IR_knp), list(topo_castle, topo_knp), f_rasterize_IR)
IR_castle <- IR_rasters[[1]]
IR_knp <- IR_rasters[[2]]


# merge -------------------------------------------------------------------

# gonna keep the two fires seperate. should lead to smaller overall space. 

# KNP ==================================

raster_list <- c(
  list(
    cbi_knp %>% setNames('cbi'),
    crop(firehistory_cbi, ndvi_knp) %>% setNames('FH_cbi'), 
    crop(firehistory_year, ndvi_knp) %>% setNames('FH_year'),
    IR_knp %>% setNames('jday'),
    ndvi_knp %>% setNames('ndvi')
  ) ,
  firDir_knp,
  topo_knp,
  crop(mtbs, ndvi_knp) %>% setNames( paste0('mtbs_', names(.)))
)

# force rasters to common extent. use ndvi layer.
rasters_knp <- map(raster_list, extend, y = ndvi_knp) %>% rast 


# Castle ==============================
raster_list <- c(
  list(
    cbi_castle %>% setNames('cbi'),
    crop(firehistory_cbi, ndvi_castle) %>% setNames('FH_cbi'), 
    crop(firehistory_year, ndvi_castle) %>% setNames('FH_year'),
    IR_castle %>% setNames('jday'),
    ndvi_castle %>% setNames('ndvi')
  ),
  firDir_castle,
  topo_castle,
  crop(mtbs, ndvi_castle) %>% setNames( paste0('mtbs_', names(.)))
)

# force rasters to common extent. use topo layer.
rasters_castle <- map(raster_list, extend, y = ndvi_castle) %>% rast 



# export ------------------------------------------------------------------

path <- 'outputs/spatial/compiled'
if(!dir.exists(path)) dir.create(path, recursive = T)
writeRaster(rasters_knp, filename = file.path(path, 'rasters_knp.tif'), overwrite = T)
writeRaster(rasters_castle, filename = file.path(path, 'rasters_castle.tif'), overwrite = T)

