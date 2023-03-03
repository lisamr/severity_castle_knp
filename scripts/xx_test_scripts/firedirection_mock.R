# how can you declare if a fire is heading, backing or flanking?
rm(list=ls())
library(sf)
library(tidyverse)
library(nngeo) #remotes::install_github("michaeldorman/nngeo")
library(terra)
theme_set(theme_classic())



#=========================================================================
# CREATE POLYGONS
#=========================================================================

# load topogrpahy
topo <- rast('outputs/topography/castle_terrain.tif')
randompt <- spatSample(topo, 1, 'random', na.rm = T, xy = T)
X0 = randompt$x
Y0 = randompt$y

# make polygons, where one is inside the other
p4 <- st_buffer(st_point(c(X0 + 1200,Y0)), 2000)
p3 <- st_buffer(st_point(c(X0 + 800,Y0)), 1000)
p2 <- st_buffer(st_point(c(X0 + 400,Y0)), 500)
p1 <- st_buffer(st_point(c(X0,Y0)), 100)
geometry <- st_sfc(p1, p2, p3, p4)
ndays <- length(geometry)
polygons <- st_sf(day = 1:ndays, geometry, crs = crs(topo)) 
#plot(arrange(polygons, rev(day)))

# crop topo
topocrop <- crop(x = topo, polygons)

#=========================================================================
# TURN POLYGONS INTO RASTER
#=========================================================================

# convert to raster
v <- vect(polygons)
pr <- rasterize(rev(v), topocrop$elevation, 'day') 


#=========================================================================
# GET AZIMUTH TO NEAREST EDGE FROM LAST DAY
#=========================================================================

# convert the pixels to points
# convert the edge of the polygons to points
# locate the nearest point from previous day
# calculate azimuth to it.

# pixels and edges to points
pr_pts <- as.points(pr) %>% st_as_sf # turn raster into points
perim_pts <- st_cast(polygons, 'LINESTRING') %>% # turn polygon edges into points
  st_line_sample(density = 1) %>% # density = points per unit (m) 
  st_sf(day = polygons$day, geometry = .) %>% 
  st_cast('POINT')


# get azimuth to nearest polygon edge
f_az_NN <- function(from, to){
  
  # is day to 1-day from? ignore for day 1
  dayfrom <- unique(from$day)
  if(dayfrom != 1 & dayfrom != unique(to$day) + 1){
    warning(glue::glue('day from is {dayfrom} and day to is {unique(to$day)}'))
  } 
  
  if(dayfrom <= 1){
    message('time = 1 will result in NA')
    return(rep(NA, nrow(from)))
  } else {
    nearest_edge <- to[st_nearest_feature(from, to),]
    azimuth <- st_azimuth(from, nearest_edge)
    return(azimuth)
  }
}

from_list <- split(pr_pts, pr_pts$day)
to_list <- split(perim_pts, perim_pts$day)
to_list <- to_list[c(1, 1:(length(to_list)-1))]
pr_pts <- map2(from_list, to_list, f_az_NN) %>% 
  map2_df(from_list, ., ~ mutate(.x, azimuth = .y, .after = day))

#=========================================================================
# DIRECTION OF FIRE
#=========================================================================

# first apply smoothing filter on aspect.
topocrop$aspect_sm <- focal(topocrop$aspect, fun = 'mean', w= 3) 

# rasterize azimuth, get direction
az_r <- rasterize(vect(pr_pts), pr, 'azimuth')
dx <- abs(az_r - topocrop$aspect) 
dx[dx > 180] <- 360 - dx
mat <- rbind(
  cbind(-Inf, 60, 1), # heading
  cbind(60, 120, 2), # flanking
  cbind(120, Inf, 3) # backing
)
dx_class <- classify(dx, mat)

fire_dir <- rast(list(az_r, dx, dx_class)) %>% 
  setNames(c('azimuth', 'dir', 'dir_class'))






# another faster way to do it ---------------------------------------------


f_direction <- function(raster, t){
  mask_from <- mask(raster, raster == t - 1, maskvalues = F)
  raster_dir <- direction(mask_from, degrees = T, from = F)
  
  mask_to <- mask(raster, raster == t, maskvalues = F)
  raster_dir <- mask(raster_dir, mask_to)
  return(raster_dir)
}

az_r2 <- lapply(2:ndays, function(t) f_direction(pr, t)) %>% sprc %>% mosaic

#par(mfrow = c(1,2))
fm <- focalMat(az_r2, 30, 'Gauss')
#plot(az_r2)
az_r2 = focal(az_r2, fm, na.rm = F)# %>% plot

# first apply smoothing filter on aspect.
topocrop$aspect_sm <- focal(topocrop$aspect, fun = 'mean', w= 3)

# rasterize azimuth, get direction
dx <- abs(az_r2 - topocrop$aspect)
dx[dx > 180] <- 360 - dx
mat <- rbind(
  cbind(-Inf, 60, 1), # heading
  cbind(60, 120, 2), # flanking
  cbind(120, Inf, 3) # backing
)
dx_class <- classify(dx, mat)

fire_dir_2 <- rast(list(az_r2, dx, dx_class)) %>%
  setNames(c('azimuth', 'dir', 'dir_class'))




#=========================================================================
# visualize
#=========================================================================

fire_dir_df <- fire_dir %>% as.data.frame(xy = T)
topo_df <- topocrop %>% as.data.frame(xy = T)

cowplot::plot_grid(
  ggplot(topo_df) +
    geom_raster(aes(x, y, fill = shade)) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_gradient(low = 'grey100', high = 'grey20'),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = azimuth)) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_gradient(low = 'grey100', high = 'grey20'),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = dir)) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_distiller(palette = 'RdYlBu', direction = 1),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = as.factor(dir_class))) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_brewer(palette = 'RdYlBu', direction = 1, 
                      labels = c('heading', 'flanking', 'backing'))
)


fire_dir_df <- fire_dir_2 %>% as.data.frame(xy = T)

cowplot::plot_grid(
  ggplot(as.data.frame(topocrop$shade, xy = T)) +
    geom_raster(aes(x, y, fill = shade)) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_gradient(low = 'grey100', high = 'grey20'),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = azimuth)) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_gradient(low = 'grey100', high = 'grey20'),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = dir)) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_distiller(palette = 'RdYlBu', direction = 1),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = as.factor(dir_class))) +
    geom_sf(data = polygons, fill = NA) +
    scale_fill_brewer(palette = 'RdYlBu', direction = 1, 
                      labels = c('heading', 'flanking', 'backing'))
)


