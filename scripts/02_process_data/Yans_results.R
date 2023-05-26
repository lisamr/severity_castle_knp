
# compile yans mortality data into a raster stack. dead canopy area should also have no NAs -> 0. 


library(terra)
library(tidyverse)


# load yan's files
yans_paths <- fs::dir_ls('data/Yan_models/', glob = '*.tif')
yans_names <- c('dead_area', 'count_grey', 'count_red', 'eccent', 'nonveg_mask')
yans_files <- map(yans_paths, ~ rast(.x))


# make all same extent/resolution
yans_files$`data/Yan_models/eccentricity_300m_EPSG32611.tif` <- 
  resample(yans_files$`data/Yan_models/eccentricity_300m_EPSG32611.tif`,
           yans_files$`data/Yan_models/area_30m_EPSG32611.tif`, method= 'bilinear')
yans_files$`data/Yan_models/mask_30m_EPSG32611.tif` <- 
  crop(yans_files$`data/Yan_models/mask_30m_EPSG32611.tif`, 
       yans_files$`data/Yan_models/area_30m_EPSG32611.tif`)
yans_files <- rast(yans_files) %>% 
  setNames(yans_names)


# convert NAs to 0 for dead area. and make it a proportion 0-1
yans_files$dead_area <- ifel(is.na(yans_files$dead_area), 0, yans_files$dead_area/900) 

# add a layer for tree counts
yans_files$dead_count <- yans_files$count_red + yans_files$count_grey


# export ------------------------------------------------------------------

path <- 'outputs/spatial/Yan'
if(!dir.exists(path)) dir.create(path, recursive = T)
writeRaster(yans_files, filename = file.path(path, 'yans_files.tif'), overwrite = T)
# 
# 
# n = 10000
# tmp <- spatSample(yans_files, n, na.rm=T)
# tmp2 <- tmp %>%
#   mutate(across(everything(), function(x){
#     (x-mean(x))/sd(x)
#   }  ))
# head(tmp2)
# v <- tmp$count_red
# HPDI(v, prob = .9); rethinking::dens(v)
# HPDI(rnorm(10000), .9)
# 
# apply(tmp, 2, function(x) mean(x ==0 ))
# 
# ggplot(tmp2, aes(count_red, dead_count)) + geom_jitter(alpha = .1)
# ggplot(tmp2, aes(count_red, count_grey)) + geom_jitter(alpha = .1)
# ggplot(tmp2, aes(dead_area, count_grey)) + geom_jitter(alpha = .1)
# ggplot(tmp2, aes(dead_area, count_red)) + geom_jitter(alpha = .1)
# ggplot(tmp2, aes(dead_area, dead_count)) + geom_jitter(alpha = .1)
# ggplot(tmp, aes(dead_area, dead_count)) + geom_jitter(alpha = .1)
# rethinking::dens(tmp2$dead_count)
