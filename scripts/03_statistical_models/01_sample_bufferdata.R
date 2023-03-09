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
N <- 1000
radii <- seq(50, 1000, by = 50) #c(30, seq(50, 1000, by = 100))
dist_dep_var <- 'ndvi'
perimeter <- fire[3,]

#============================
#============================


# assign sample locations -------------------------------------------------


# sample locations
set.seed(1)
loc <- spatSample(mask(r_brick$ndvi, perimeter), N, 
                  xy = T, as.df = T, na.rm = T) %>% 
  as_tibble() %>% 
  mutate(siteID = 1:N)



# make buffers instead of rings -------------------------------------------

buffers <- f_buffers(loc, radii) 

# extract mean and sd from raster
f <- function(x, na.rm = T) {
  c(mean=mean(x, na.rm = na.rm),
    sd=sd(x, na.rm = na.rm)
  )
}

# split data into 16 chunks, pass to furrr
ncores <- parallel::detectCores()
chunks <- rep(1:ncores, length.out = nrow(buffers)) %>% sort
buffers$chunks <- chunks
buffers_split <- split(buffers, chunks) %>% map(., ~ mutate(.x, inchunk_ID = 1:nrow(.x)))
buffers <- bind_rows(buffers_split)

# do in foreach?
cl <- parallel::makeCluster(ncores - 1)
doParallel::registerDoParallel(cl)
tic()
ndvi_packed <- wrap(r_brick[dist_dep_var]) # terra requires packing and unpacking in parrallel
res <- foreach(i = 1:max(chunks)) %dopar% {
  r <- terra::unwrap(ndvi_packed)
  terra::extract(r, terra::vect(buffers_split[[i]]), fun = f)
}
toc()
parallel::stopCluster(cl)

# merge back together
res <- map(res, ~ rename(.x, inchunk_ID = ID))
metrics_b <- map2(buffers_split, res, full_join)  %>% 
  bind_rows() %>% 
  select(-c(chunks, inchunk_ID)) 
df_buffers <- metrics_b %>% 
  st_drop_geometry() 
  

# sample rest of rasters --------------------------------------------------

preds <- terra::extract(r_brick, loc[,c('x', 'y')])

dir_path <- 'outputs/tabular'
if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
write_rds(list(buffers = df_buffers, preds = preds), file.path(dir_path, 'model_knp_bufferdata.rds'))
