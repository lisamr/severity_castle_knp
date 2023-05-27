# lay grid over castle and knp fires

library(tidyverse)
library(sf)


# =================
# set arguments
# =================

grid_width <- 180


# fire perimeters
fires <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')
knp <- fires[3,]
castle <- fires[1,]


# define sampling grid 
grid_knp <- st_make_grid(knp, what = 'centers', cellsize = grid_width, crs = st_crs(fires))
grid_knp <- grid_knp[st_intersects(grid_knp, knp, sparse = F)]

grid_castle <- st_make_grid(castle, what = 'centers', cellsize = grid_width, crs = st_crs(fires))
grid_castle <- grid_castle[st_intersects(grid_castle, castle, sparse = F)]


# export it
if(!dir.exists('outputs/spatial/compiled/')) dir.create('outputs/spatial/compiled/')
write_rds(grid_knp, 'outputs/spatial/compiled/grid_knp.rds')
write_rds(grid_castle, 'outputs/spatial/compiled/grid_castle.rds')
