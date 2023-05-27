# for each focal location i, calculate mean value of raster X in non-overlapping
# concentric rings of radius r
rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)

source('scripts/functions/functions_process_spatial.R')

# fire footprint
fire <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')


#============================
# assign parameters
#============================
r_brick <- rast('outputs/spatial/compiled/rasters_knp.tif') 
N <- 1000
radii <- seq(50, 2000, by = 100) #c(30, seq(50, 1000, by = 100))
dist_dep_var <- 'ndvi'
perimeter <- fire[3,]

#============================
#============================


# assign sample locations -------------------------------------------------


# sample locations
loc <- spatSample(mask(r_brick$ndvi, perimeter), N, 
                  xy = T, as.df = T, na.rm = T) %>% 
  as_tibble() %>% 
  mutate(id = 1:N)


# distance dependent variable ---------------------------------------------


# create concentric rings around location i
rings <- f_rings(loc, radii, st_crs(perimeter))

# extract mean and sd of chosen raster in each ring
f <- function(x, na.rm = T) {
  c(mean=mean(x, na.rm = na.rm),
    sd=sd(x, na.rm = na.rm)
  )
}
metrics <- terra::extract(r_brick[dist_dep_var], vect(rings), fun = f) 
metrics <- metrics %>% 
  select(-ID) %>% 
  set_names(c('mean', 'sd'))

df_rings <- rings %>% 
  mutate(ha = st_area(.) %>% units::set_units(ha)) %>% 
  as_tibble() %>% 
  select(-geometry) %>% 
  mutate(metrics) 
  
  

# sample rest of rasters --------------------------------------------------

preds <- terra::extract(r_brick, loc[,c('x', 'y')])

dir_path <- 'outputs/tabular'
if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
write_rds(list(rings = df_rings, preds = preds), file.path(dir_path, 'model_knp_data.rds'))
