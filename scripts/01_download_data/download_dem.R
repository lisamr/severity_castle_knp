# download dem for knp and castle fires

# load packages
library(sf)
library(elevatr)
library(tidyverse)
library(furrr)
library(fs)
library(raster)

# load fires
perims <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')

perims_list <- split(perims, 1:nrow(perims)) %>% 
  map(st_buffer, 1500) %>% # buffer by more than the largest window of TPI
  map(st_transform, 4326) # unprojected needed for getting dem

# use future bc this takes a bit of time.
plan(multisession, workers = length(perims_list))
dem_list <- future_map(perims_list, 
           ~ get_elev_raster(locations = .x, clip = 'locations', z = 13) %>% 
             setNames(., 'elevation')
           )
plan(sequential)
names(dem_list) <- perims$FIRE_NAME

dir_path <- 'outputs/spatial/topography'
if(!dir_exists(dir_path)) dir_create(dir_path)
write_rds(dem_list, file.path(dir_path, 'dem_list.rds'))

