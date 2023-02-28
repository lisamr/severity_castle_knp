# derive metrics from dem, put into raster brick:
# elevation, slope, aspect, tpi, hillshade, hli, concavity, etc.

library(terra)
library(tidyverse)
library(raster) # needed for loading dem
library(spatialEco)
library(glue)

# use fire severity data as template for projection
cbi <- rast('outputs/GEE_CBI/fireHistory_1985_2020_CBI_bc.tif') 

# load dem data, convert to terra spatRast
dem_list <- read_rds('outputs/topography/dem_list.rds') 
dem_list <- map(dem_list, rast)

# reproject rasters
dem_list_proj <- map(dem_list, 
                     ~ project(.x, cbi, method = 'bilinear', align = T) %>% 
                       setNames('elevation'))



# get terrain metrics -----------------------------------------------------


f_terrain_metrics <- function(elev){
  # slope, aspect
  slope_aspect <- terrain(elev, c('slope', 'aspect'))
  
  # HLI
  heatload <- hli(elev) %>% setNames('HLI')
  
  # TPI at various windows
  dist_vector <- c(100, 250, 500, 1000)
  TPI_list <- map(dist_vector, 
                  ~ tpi(elev, win = 'circle', scale = .x, normalize = T) %>% 
                    setNames(., glue('TPI_{.x}')))
  
  
  # join into a raster brick
  brick <- rast(c(list(elev, slope_aspect, heatload), TPI_list))
  return(brick)
  
}

terrain_list <- map(dem_list_proj, f_terrain_metrics)



# export ------------------------------------------------------------------


fire_names <- sapply(strsplit(names(terrain_list), " "), `[`, 1) %>% 
  janitor::make_clean_names()

dir_path <- 'outputs/topography'
if(!fs::dir_exists(dir_path)) fs::dir_create(dir_path)
walk(seq_along(fire_names), 
     ~ writeRaster(terrain_list[[.x]], glue('{dir_path}/{fire_names[.x]}_terrain.tif')))




