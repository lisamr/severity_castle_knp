# create sample locations. KNP only for now.
# sample from spatiotemporal data (fire dates, climate)

rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)
library(foreach)

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_sample_bricks.R')

# grid locations
grid_loc <- read_rds('outputs/spatial/compiled/grid_knp.rds') 

# perimeters from IR with VIIRS filled in
perimeters <- read_rds('outputs/spatial/IR/knp_growth_IR_VIIRS.rds')


# weather data
nc_paths <- fs::dir_ls('data/Williams_weather/', recurse = T, glob = '*.nc')
nc_names <- sub("\\.nc$", "", basename(nc_paths))
wx_list <- map(nc_paths, ~ rast(.x)) %>% 
  set_names(nc_names)

# raster to resample to
r <- rast('outputs/spatial/Yan/yans_files.tif')


# 
# # sample climate ----------------------------------------------------------
# 
# # climate data too (sample with points rather than reproject)
# climate <- rast('outputs/spatial/weather/climate_means_2000_2020_MayOct.tif')
# 
# # convert points to crs of climate
# grid_loc_lonlat <- grid_loc %>% 
#   st_transform(crs = 'EPSG:4326') %>% 
#   st_sf
# 
# # sample climate
# climate_extract <- terra::extract(climate, grid_loc_lonlat)
# climate_pts <- grid_loc %>% 
#   st_sf() %>% 
#   mutate(climate_extract, .before = geometry)
# 
# grid_loc %>% 
#   st_sf() %>% 
#   ggplot(.) +
#   geom_sf()
# 
# ggplot(climate_pts, aes(color = rh)) + geom_sf()

# get burn dates ----------------------------------------------------------

# extract burn date ~ 20sec
grid_locsf <- st_sf(grid_loc)
perimeters <- st_transform(perimeters, st_crs(grid_locsf)) 

int_matrix <- st_intersects(grid_locsf, perimeters, sparse = F) 
int_idx <- apply(int_matrix, 1, function(x) which(x)[1])
burn_dates <- as_date(perimeters$date_adj[int_idx])

grid_loc_filt <- grid_locsf %>% 
  mutate(date_adj = burn_dates) %>% 
  filter(!is.na(date_adj))


# extract weather ---------------------------------------------------------


# shorten the climate data to only the burn dates. runs faster. 
sampled_dates <- grid_loc_filt$date_adj %>% unique %>% sort

system.time(
  wx_sampled <- map(wx_list, ~ subset_brick(.x, sampled_dates)) %>% 
    map(~ project(.x, r)) %>% 
    map_df(~ sample_brick(grid_loc_filt, .x)) 
) # ~10s without reprojecting, 144 with

wx_sampled_sf <- bind_cols(grid_loc_filt, wx_sampled)

# check to make sure they make sense
ggplot(wx_sampled_sf, aes(color = wind)) +
  geom_sf()
ggplot(grid_loc_filt, aes(color = date_adj)) +
  geom_sf()

# export ------------------------------------------------------------------

write_rds(wx_sampled_sf, 'outputs/spatial/weather/spatiotemporal_weather_KNP.rds')



# visualize ---------------------------------------------------------------


# temp/vpd very similar. rh, wind, solar good candidates. 
wx_sampled_sf %>% 
  st_drop_geometry() %>% 
  sample_n(1000) %>% pairs

# rh/vpd and wind seem like important vars
ggplot(wx_sampled_sf, aes(color = wind)) +
  geom_sf()
ggplot(wx_sampled_sf, aes(color = vpd)) +
  geom_sf()
ggplot(wx_sampled_sf, aes(color = solar)) +
  geom_sf()
ggplot(wx_sampled_sf, aes(color = rh)) +
  geom_sf()
ggplot(wx_sampled_sf, aes(color = tmax)) +
  geom_sf()
ggplot(wx_sampled_sf, aes(color = tmin)) +
  geom_sf()
ggplot(wx_sampled_sf, aes(color = prec)) +
  geom_sf()


