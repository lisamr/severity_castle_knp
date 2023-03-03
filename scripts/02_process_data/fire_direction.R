# how can you declare if a fire is heading, backing or flanking?
rm(list=ls())
library(sf)
library(tidyverse)
library(nngeo) #remotes::install_github("michaeldorman/nngeo")
library(terra)
library(tictoc)
theme_set(theme_classic())



# load data ---------------------------------------------------------------


# load topogrpahy
topo_castle <- rast('outputs/spatial/topography/castle_terrain.tif')
topo_knp <- rast('outputs/spatial/topography/knp_terrain.tif')

# load IR
IR_castle <- st_read('outputs/spatial/IR/castle_byDay_u.geojson')
IR_knp <- st_read('outputs/spatial/IR/knp_byDay_u.geojson')




# get fire direction ------------------------------------------------------


f_fire_direction <- function(IR_data, topo_data){
  # check data inputs
  stopifnot(class(topo_data)[1] == 'SpatRaster')
  stopifnot(st_geometry_type(IR_data) %>% unique %in% c('MULTIPOLYGON', 'POLYGON'))
  
  # rasterize IR data -------------------------------------------------------
  
  # add consecutive day index
  IR_data$dayID <- 1:nrow(IR_data)
  
  # rasterize IR data
  IR_data_r <- rasterize(rev(vect(IR_data)), topo_data$elevation, 'dayID') 
  
  
  # azimuth to nearest edge  ------------------------------------------------
  
  # pixels and edges to points
  IR_data_pts <- as.points(IR_data_r) %>% st_as_sf # turn raster into points
  tic()
  casted <- st_cast(IR_data, 'POLYGON') %>% # turn polygon edges into points
    st_cast('LINESTRING') 
  message('IR data casted')
  perim_pts <- casted %>% 
    st_line_sample(density = 1/10) %>% 
    st_sf(dayID = casted$dayID, geometry = .) %>% 
    st_cast('POINT') %>% 
    filter(!st_is_empty(.))
  message('IR data turned into points')
  toc()
  
  
  # get azimuth to nearest polygon edge
  f_az_NN <- function(from, to){
    tic()
    # is day to 1-day from? ignore for day 1
    dayfrom <- unique(from$dayID)
    if(dayfrom != 1 & dayfrom != unique(to$dayID) + 1){
      warning(glue::glue('day from is {dayfrom} and day to is {unique(to$dayID)}'))
    } 
    
    if(dayfrom <= 1){
      message('time = 1 will result in NA')
      return(rep(NA, nrow(from)))
    } else {
      nearest_edge <- to[st_nearest_feature(from, to),]
      azimuth <- st_azimuth(from, nearest_edge)
      cat('azimuth finished for dayID: ', dayfrom, '\n')
      print(toc())
      return(azimuth)
    }
  }
  
  from_list <- split(IR_data_pts, IR_data_pts$dayID)
  to_list <- split(perim_pts, perim_pts$dayID)
  # ensure that there's exactly 1 day between to and from
  dayIDs_from <- from_list %>% map_vec(~unique(.x$dayID))
  dayIDs_to <- to_list %>% map_vec(~unique(.x$dayID))
  to_list <- to_list[dayIDs_to[c(1, dayIDs_from-1)]]
  IR_data_pts <- map2(from_list, to_list, f_az_NN) %>% 
    map2_df(from_list, ., ~ mutate(.x, azimuth = .y, .after = dayID))
  
  # direction of fire -------------------------------------------------------
  
  
  # first apply smoothing filter on aspect.
  topo_data$aspect_sm <- focal(topo_data$aspect, fun = 'mean', w= 3) 
  
  # rasterize azimuth, get direction
  az_r <- rasterize(vect(IR_data_pts), IR_data_r, 'azimuth')
  dx <- abs(az_r - topo_data$aspect_sm) 
  dx[dx > 180] <- 360 - dx
  mat <- rbind(
    cbind(-Inf, 60, 1), # heading
    cbind(60, 120, 2), # flanking
    cbind(120, Inf, 3) # backing
  )
  dx_class <- classify(dx, mat)
  
  fire_dir <- rast(list(az_r, dx, dx_class)) %>% 
    setNames(c('azimuth', 'dir', 'dir_class'))
  
  return(fire_dir)
}

fire_dir_knp <- f_fire_direction(IR_knp, topo_knp)
fire_dir_castle <- f_fire_direction(IR_castle, topo_castle)



#=========================================================================
# EXPORT
#=========================================================================

dir_path <- 'outputs/spatial/fire_direction'
if(!dir.exists(dir_path)) dir.create(dir_path)
writeRaster(fire_dir_castle, file.path(dir_path, 'fire_dir_castle.tif'))
writeRaster(fire_dir_knp, file.path(dir_path, 'fire_dir_knp.tif'))


#=========================================================================
# visualize
#=========================================================================

fire_dir <- fire_dir_knp
fire_dir_df <- fire_dir %>% as.data.frame(xy = T)
shade_df <- as.data.frame(topo_knp$shade, xy = T)
IR_data <- IR_knp


# see all 
cowplot::plot_grid(
  ggplot(shade_df) +
    geom_raster(aes(x, y, fill = shade)) +
    geom_sf(data = IR_data, fill = NA) +
    scale_fill_gradient(low = 'grey100', high = 'grey20'),
  ggplot(fire_dir_df) +
    geom_raster(aes(x, y, fill = dir)) +
    geom_sf(data = IR_data, fill = NA) +
    scale_fill_distiller(palette = 'RdYlBu', direction = 1)
)


# zoom in
rand_pt <- spatSample(fire_dir, 1, xy = T, na.rm = T)[,c('x', 'y')]
xlims <- c(rand_pt$x - 2000, rand_pt$x + 2000)
ylims <- c(rand_pt$y - 2000, rand_pt$y + 2000)
shade_df_f <- filter(shade_df, 
                     x >= xlims[1] & x <= xlims[2] & 
                     y >= ylims[1] & y <= ylims[2]  )
fire_dir_df_f <- filter(fire_dir_df, 
                     x >= xlims[1] & x <= xlims[2] & 
                       y >= ylims[1] & y <= ylims[2]  )

cowplot::plot_grid(
  # shade
  ggplot(shade_df_f) +
    geom_raster(aes(x, y, fill = shade)) +
    geom_sf(data = IR_data, fill = NA, aes(color = jday)) +
    scale_fill_gradient(low = 'grey100', high = 'grey20') +
    scale_color_viridis_c() +
    coord_sf(crs = st_crs(IR_data), xlim = xlims, ylim = ylims) +
    theme(legend.position = 'none') ,
  # direction
  ggplot(fire_dir_df_f) +
    geom_raster(aes(x, y, fill = dir)) +
    geom_sf(data = IR_data, fill = NA, aes(color = jday)) +
    scale_fill_distiller(palette = 'RdYlBu', direction = 1)+
    scale_color_viridis_c() +
    coord_sf(crs = st_crs(IR_data), xlim = xlims, ylim = ylims) +
    theme(legend.position = 'none'),
  # direction class
  ggplot(fire_dir_df_f) +
    geom_raster(aes(x, y, fill = as.factor(dir_class))) +
    geom_sf(data = IR_data, fill = NA, aes(color = jday)) +
    scale_color_viridis_c() +
    scale_fill_brewer(palette = 'RdYlBu', direction = 1, 
                      labels = c('heading', 'flanking', 'backing')) +
    coord_sf(crs = st_crs(IR_data), xlim = xlims, ylim = ylims) +
    theme(legend.position = 'none'),
  nrow = 1
)



# trash -------------------------------------------------------------------



# another way to do it--use `terra::direction`  ---------------------------
# 
# 
# f_direction <- function(raster, tm){
#   tic()
#   mask_from <- mask(raster, raster == tm - 1, maskvalues = F)
#   raster_dir <- direction(mask_from, degrees = T, from = F)
#   
#   mask_to <- mask(raster, raster == tm, maskvalues = F)
#   raster_dir <- mask(raster_dir, mask_to)
#   cat('finished t:', tm, '\n')
#   print(toc())
#   return(raster_dir)
# }
# dayIDs <- sort(unique(IR_data$dayID))
# fire_az <- lapply(dayIDs[-1], function(tm) f_direction(IR_data_r, tm)) %>% 
#   sprc %>% mosaic
# 
# # first apply smoothing filter on aspect.
# topo_data$aspect_sm <- focal(topo_data$aspect, fun = 'mean', w= 3)
# 
# # rasterize azimuth, get direction
# dx <- abs(fire_az - topo_data$aspect_sm)
# dx[dx > 180] <- 360 - dx
# mat <- rbind(
#   cbind(-Inf, 60, 1), # heading
#   cbind(60, 120, 2), # flanking
#   cbind(120, Inf, 3) # backing
# )
# dx_class <- classify(dx, mat)
# 
# fire_dir <- rast(list(fire_az, dx, dx_class)) %>%
#   setNames(c('azimuth', 'dir', 'dir_class'))
# 
# 
# 




