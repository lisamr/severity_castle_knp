#install.packages('ncdf4')

library(terra)
library(tidyverse)
library(fs)
library(sf)
library(lubridate)

# can read in .nc files with terra
nc_paths <- dir_ls('data/Williams_weather/', recurse = T, glob = '*.nc')
nc_names <- sub("\\.nc$", "", basename(nc_paths))


# how would you extract values at a certain location on a certain date? 
knp <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp') 
knp <- knp[,3]

# make an sf dataset
t0 <- as_date('2021-09-10') 
knp_dates <- t0 + 0:9
n = 25
pts <- st_sample(knp, n) %>% 
  st_sf() %>% 
  mutate(date = sample(knp_dates, n, T), .before = geometry)

# sample raster at those points at that time
sample_brick <- function(pts, brick){
  pts_unproj <- pts %>% arrange(date) %>% st_transform(crs(brick))
  date_var <- time(brick)
  col_idx <- match(pts_unproj$date, date_var)
  
  # extract all data. rows=locations, cols = dates
  s <- terra::extract(brick, pts_unproj, ID = F, layer = col_idx)

  return(s$value)
}

tic()
sample_brick(pts, tmp)
toc() # 6.64 sec elapsed

smaller_tmp <- tmp[[which(time(tmp) %in% knp_dates)]]
tic()
s <- sample_brick(pts, smaller_tmp)
toc() # 6.64 sec elapsed
s$value



