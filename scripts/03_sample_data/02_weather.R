# create sample locations. 
# sample from spatiotemporal data (fire dates, climate)

rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)
library(foreach)

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_sample_bricks.R')

# fire footprint
fire <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')
knp <- fire[3,]

# VIIRS-derived perimeters
viirs_list <- read_rds('outputs/spatial/VIIRS/perims_viirs.rds')

# weather data
nc_paths <- fs::dir_ls('data/Williams_weather/', recurse = T, glob = '*.nc')
nc_names <- sub("\\.nc$", "", basename(nc_paths))
wx_list <- map(nc_paths, ~ rast(.x)) %>% 
  set_names(nc_names)
 
# assign sample locations -------------------------------------------------

N <- 100
set.seed(1)
pts_knp <- st_sample(knp, N) %>% st_sf()

# extract burn date 
burn_dates <- suppressWarnings(
  st_intersection(pts_knp, viirs_list$polys_growth)
) 


# extract from climate data -----------------------------------------------

# shorten the climate data to only the burn dates. runs faster. 
sampled_dates <- burn_dates$date %>% unique
wx_sampled <- map(wx_list, ~ subset_brick(.x, sampled_dates)) %>% 
  map_df(~ sample_brick(burn_dates, .x)) %>% 
  mutate(ID = 1:nrow(.), date = burn_dates$date, .before = prec)

wx_sampled_sf <- st_sf(wx_sampled, geometry = st_geometry(burn_dates))


