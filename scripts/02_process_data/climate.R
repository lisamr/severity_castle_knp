# get long term climate data

rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)
library(foreach)
library(lubridate)
library(parallel)

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_sample_bricks.R')

# weather data
nc_paths <- fs::dir_ls('data/Williams_weather/', recurse = T, glob = '*.nc')
nc_names <- sub("\\.nc$", "", basename(nc_paths))
wx_list <- map(nc_paths, ~ rast(.x)) %>% 
  set_names(nc_names)

# dates to consider
startdate <- as_date('2000-05-01') 
enddate <- as_date('2020-09-01') 
date_sequence <- seq(startdate, enddate, by = 'day')
# only consider fire season (may1 - oct31)
date_sequence <- date_sequence[month(date_sequence) >= 5 & month(date_sequence) <= 10]

# subset the data to those dates
wx_list_subset <- map(wx_list, ~ subset_brick(.x, date_sequence))


# get averages for every pixel and climate variable. 
# takes a while. can't figure out how to parallelize terra objects with e.g. foreach
wx_mean <- vector('list', length(wx_list_subset))
t1 <- proc.time()
for(i in 1:length(wx_list_subset)){
  print(paste0('starting i=', i))
  tic()
  wx_mean[[i]] <- app(wx_list_subset[[i]], mean, cores = 15)
  toc()
  print(paste0('finished i=', i))
}
t2 <- proc.time()
ttotal <- t2-t1
ttotal[3]/60 # 67.4 minutes

# merge into a single brick
climate_brick <- rast(wx_mean)
names(climate_brick) <- nc_names

# export
writeRaster(climate_brick, 'outputs/spatial/weather/climate_means_2000_2020_MayOct.tif')
