# look at the files from michele. should be live trees from 2019, modeled from F3.


rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)

# look at available files
fs::dir_ls('data/F3_2019_livetrees/', glob = '*.tif') %>% basename()
# [1] "Total3Run_BA_24UP_NoMGT_2019_V20220512.tif" # live basal area, 24in+ dbh, no mgt F3 trt
# [2] "Total3Run_TPA24UP_NoMGT_2019_V20220512.tif" # live tree counts, 24in+ dbh, no mgt F3 trt
# [3] "Total3Run_TPA_NoMGT_2019_V20220512.tif" # live tree counts, no mgt F3 trt


# needs to be consistent with the rest of the rasters. use yans data.
yans_files <- rast('outputs/spatial/Yan/yans_files.tif')


TPA <- rast('data/F3_2019_livetrees/Total3Run_TPA24UP_NoMGT_2019_V20220512.tif')
TPA_reproj <- project(TPA, yans_files)

# export
if(!dir.exists('outputs/spatial/livetrees')) dir.create('outputs/spatial/livetrees')
writeRaster(TPA_reproj, 'outputs/spatial/livetrees/TPA24.tif')
