# for each focal location i, calculate mean value of raster X in non-overlapping
# concentric rings of radius r
rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)
library(foreach)

source('scripts/functions/functions_process_spatial.R')

# fire footprint
fire <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')


#============================
# assign parameters
#============================
r_brick <- rast('outputs/spatial/compiled/rasters_knp.tif') 
weather <- read_rds('outputs/spatial/weather/spatiotemporal_weather_KNP.rds')

radii <- seq(50, 1000, by = 50) #c(30, seq(50, 1000, by = 100))
dist_dep_var <- c('dead_area', 'count_red', 'count_grey', 'TPA')
perimeter <- fire[3,]


#============================
#============================

# sample locations
loc <- st_coordinates(weather)
colnames(loc) <- c('x', 'y')





# make rings --------------------------------------------------------------

# user  system elapsed 
# 231.53   12.14  244.40 
system.time(rings <- f_rings(loc, radii, st_crs(fire)) )

tic()
df_rings <- bind_cols(rings, 
          terra::extract(r, vect(rings), ID = F, fun = function(x) mean(x, na.rm = T))
          ) 
toc() #233.68 sec elapsed


# sample rest of rasters --------------------------------------------------

# predictors from raster
preds <- terra::extract(r_brick, loc) %>% 
  rename_with(~ paste0(., '_avg'), prec:wind) 

# merge with weather predictors
all_vars <- weather %>% 
  rename_with(~ paste0(., '_day'), prec:wind) %>% 
  mutate(preds)


dir_path <- 'outputs/tabular'
if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
write_rds(list(loc = loc, rings = df_rings, all_vars = all_vars), file.path(dir_path, 'model_knp_ringdata.rds'))




