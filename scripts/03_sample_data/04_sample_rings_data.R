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
names(r_brick)
dist_dep_var <- c('dead_area', 'dead_count', 'count_red', 'count_grey', 'TPA', 'cover')
dist_dep_var_conifer <- c('dead_area_con', 'dead_count_con', 'count_red_con', 
                          'count_grey_con', 'TPA_conifer', 'cover_conifer')

perimeter <- fire[3,]


#============================
#============================

# sample locations
loc <- st_coordinates(weather)
colnames(loc) <- c('x', 'y')





# make rings --------------------------------------------------------------

r <- r_brick[[dist_dep_var]]
r_conifer <- r_brick[[dist_dep_var_conifer]]

# user  system elapsed 
# 231.53   12.14  244.40 
system.time(rings <- f_rings(loc, radii, st_crs(fire)) )

#rings <- read_rds('outputs/tabular/rings.rds')

# extract with exact_extract. orders of magnitude faster than terra::extract
tic()
rings_ext <- exactextractr::exact_extract(x = r, y = rings, fun = "mean")
toc() #106.16 sec elapsed!!!!!
rings_ext_conifer <- exactextractr::exact_extract(x = r_conifer, y = rings, fun = "mean")

# rename the columns
names(rings_ext) <- str_replace(names(rings_ext), "^mean\\.", "")
names(rings_ext_conifer) <- glue::glue("{names(rings_ext)}_con")

# bind with the ring polygons
rings_df <- bind_cols(rings, rings_ext)
rings_df_conifer <- bind_cols(rings, rings_ext_conifer)

# sample rest of rasters --------------------------------------------------

# predictors from raster
preds <- terra::extract(r_brick, loc) %>% 
  rename_with(~ paste0(., '_avg'), prec:wind) 

# merge with weather predictors
all_vars <- weather %>% 
  rename_with(~ paste0(., '_day'), prec:wind) %>% 
  mutate(preds)

all_vars_df <- st_drop_geometry(all_vars) %>% 
  mutate(geometry = st_geometry(all_vars))

dir_path <- 'outputs/tabular'
if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
write_rds(list(loc = loc, rings = rings_df, rings_conifer = rings_df_conifer,
               all_vars = all_vars_df), 
          file.path(dir_path, 'model_knp_ringdata.rds'))
#write_rds(rings, file.path(dir_path, 'rings.rds'))



# # extract data in parallel ------------------------------------------------
# 
# 
# rings_split <- rings %>% 
#   mutate(chunk = sort(rep(1:16, length.out = nrow(.)))) %>% 
#   group_split(chunk) %>% 
#   map(~ mutate(.x, withinchunk_id = 1:nrow(.x)))
# 
# # convert terra r into raster r
# library(raster)
# r_raster <- brick(r)
# 
# library(furrr)
# plan(multisession)
# 
# small_rings_split <- rings_split %>% map(~sample_n(.x, 1)) %>% 
#   map(~ as(.x, 'Spatial'))
# 
# 
# my_function <- function(r, v) {
#   raster::extract(r, v, fun = mean)
# }
# library(parallel)
# library(doParallel)
# cl <- makeCluster(15)
# registerDoParallel(cl)
# out=clusterEvalQ(cl,library(raster))
# env= where('my_function')
# clusterExport(cl, list('my_function'))
# 
# system.time({
#   res <- foreach(l = c(1:16)#,.errorhandling = 'remove'
#                  ) %do%
#     my_function(r_raster, small_rings_split[[l]])
# }) 
# 
