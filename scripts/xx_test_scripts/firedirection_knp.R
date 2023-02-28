# apply fire direction script to real fire perimeters
# things to play with: density of points along the fire perimeters, smoothing of aspect layer
rm(list=ls())
library(sf)
library(tidyverse)
library(fasterize)
library(nngeo) #remotes::install_github("michaeldorman/nngeo")
library(cowplot)
library(raster)
library(lubridate)
theme_set(theme_classic())


#=========================================================================
# IMPORT AND CLEAN 
#=========================================================================

# import fire perimeters
perimeters <- st_read('data/KNP_0914_1016.shp')
#plot(perimeters$geometry)

# load topography
topo <- brick('../../../Documents/ArcGIS/Spatial_data/from_GEE/topography30m.tif')

# inspect dates of the perimeters
as_date(perimeters$CreateDate) %>% unique %>% sort


# clean up the data
perimeters <- perimeters %>% 
  mutate(# make column values consistent
         IncidentNa = 'KNPComplex',
         # add a variable for day order
         day = as.numeric(as.factor(date(CreateDate)))
         ) 
perimeters %>% 
  as.data.frame() %>% 
  distinct(CreateDate, day) %>% 
  arrange(day)

# transform to raster crs
perimeters_utm <- st_transform(perimeters, crs(topo))

# crop topography
topocrop <- crop(topo, perimeters_utm)

plot(topocrop$hillshade, col = grey.colors(10))
plot(perimeters_utm$geometry, add = T, border = sf.colors())



#=========================================================================
# TURN POLYGONS INTO RASTER
#=========================================================================

# convert to raster
perimeters_r <- fasterize(perimeters_utm, topocrop$elevation, field = 'day', fun = 'min')
names(perimeters_r) <- 'day'

plot(topocrop$hillshade, col = grey.colors(10), legend = F)
plot(perimeters_r, add = T, col = hcl.colors(20))
plot(perimeters_utm$geometry, add = T)

# export raster for Danny
writeRaster(perimeters_r, 'outputs/KNP_dailyperim.tif')

#=========================================================================
# GET AZIMUTH TO NEAREST EDGE FROM LAST DAY
#=========================================================================

# how do you calculate azimuth to nearest pixel of given value?
#https://stackoverflow.com/questions/60990366/calculating-azimuth-of-nearest-features

# convert the pixels to points
# convert the edge of the polygons to points
# locate the nearest point from previous day
# calculate azimuth to it.

# pixels and edges to points
perimeters_pts <- rasterToPoints(perimeters_r, spatial = T) %>% st_as_sf() %>% arrange(day) #raster as points
edges_lines <- st_cast(perimeters_utm, 'POLYGON') %>% st_cast('LINESTRING')  #edges as lines, ignore warning?
edges_pts <- st_line_sample(edges_lines, density= 1/500) %>% #convert to multipoints
  st_sf(day = edges_lines$day, .) %>% #add attributes
  st_cast('POINT') #convert to points
edges_pts$ID <- 1:nrow(edges_pts) #add ID to perimeter points for keeping track

# confirm that days of edges are named correctly
tmp = edges_pts[sample(nrow(edges_pts), 10000),]
ggplot() +
  geom_sf(data = tmp, aes(color = day)) +
  scale_color_viridis_c()
tmp = perimeters_pts[sample(nrow(perimeters_pts), 10000),]
ggplot() +
  geom_sf(data = tmp, aes(color = day)) +
  scale_color_viridis_c()


#find the nearest edge and add to sf object
nearestedge <- rep(NA, nrow(perimeters_pts))
v1 <- 1
ndays <- max(perimeters$day)
for(i in 1:ndays){
  index <- st_nearest_feature(filter(perimeters_pts, day == i), filter(edges_pts, day == i-1))
  ID <- filter(edges_pts, day == i-1)[index,] %>% pull('ID')
  v2 <- length(index) -1 + v1
  nearestedge[v1:v2] <- ID 
  v1 <- v2 + 1
}
perimeters_pts$nearestedge <- nearestedge
nearestedge1 <- nearestedge
nearestedge1[is.na(nearestedge1)] <- 1 #st_azimuth doesn't allow NAs, so temporarily assign NA values as 1
perimeters_pts$nearestedge1 <- nearestedge1


#calculate azimuth (and distance so I know I'm doing it right) to nearest edge
perimeters_pts$dist <- st_distance(perimeters_pts, edges_pts[nearestedge,], by_element = T) %>% as.numeric()
perimeters_pts$az <- st_azimuth(perimeters_pts, edges_pts[nearestedge1,])
perimeters_pts$az[is.na(nearestedge)] <- NA #change the values within first polygon as NA


# viz. I think that looks right.
perimeter_sample <- perimeters_pts[seq(1, nrow(perimeters_pts), length = 10000), ]
plot1 <- ggplot() +
  geom_sf(data = perimeter_sample, aes(color = az), size = 1) +
  geom_sf(data = st_geometry(perimeters_utm), size = 1, fill = NA, color = 'black') +
  scale_color_viridis_c(na.value = NA)
plot2 <- ggplot() +
  geom_sf(data = perimeter_sample, aes(color = log(dist)), size = 1) +
  geom_sf(data = st_geometry(perimeters_utm), size = 1, fill = NA, color = 'black') +
  scale_color_viridis_c(na.value = NA, )
plot3 <- ggplot() +
  geom_sf(data = perimeter_sample, aes(color = nearestedge), size = 1) +
  geom_sf(data = st_geometry(perimeters_utm), size = 1, fill = NA, color = 'black') +
  scale_color_viridis_c(na.value = NA)
plot4 <- ggplot() +
  geom_sf(data = perimeter_sample, aes(color = day), size = 1) +
  geom_sf(data = st_geometry(perimeters_utm), size = 1, fill = NA, color = 'black') +
  scale_color_viridis_c(na.value = NA)
plot_grid(plot1, plot2, plot3, plot4)




#=========================================================================
# DIRECTION OF FIRE
#=========================================================================

# first apply smoothing filter on aspect.
topocrop$aspect_sm <- focal(topocrop$aspect, w=matrix(1/9,nrow=3,ncol=3)) 

# get topo values for all the raster points
topovals <- extract(topocrop, perimeters_pts)
perimeters_pts2 <- cbind(perimeters_pts, topovals)

# get fire direction
perimeters_pts2 <- perimeters_pts2 %>%
  mutate(DX = abs(az - aspect_sm),
         DX = ifelse(DX>180, 360-DX, DX),
         DX_class = case_when( #needs to be numerical to be converted to raster
           DX <= 60 ~ 3, #heading
           DX <= 120 ~ 2, #flanking
           T ~ 1 #backing
         ),
         DX_class = ifelse(is.na(DX), NA, DX_class))

# convert to raster layer, get real dates
library(stars)
perimeters_r2 <- st_rasterize(perimeters_pts2, st_as_stars(perimeters_r))
date_tab <- as.data.frame(perimeters_utm) %>% 
  distinct(day, CreateDate)
date_tabv <- date_tab$CreateDate
names(date_tabv) <- date_tab$day
perimeters_r2 <- perimeters_r2 %>% 
  mutate(date = recode(day, !!!date_tabv))

# visualize
topocropdf <- as.data.frame(topocrop, xy = T)
plot5 <- ggplot() +
  geom_stars(data = perimeters_r2['DX']) +
  #geom_sf(data = perimeters_utm, fill = NA, color = grey(.5), lwd = .5) +
  scale_fill_viridis_c(na.value = NA, direction = -1) +
  geom_sf(data = perimeters_utm, fill = NA, aes(color = day), lwd = .2, show.legend = F) +
  scale_color_viridis_c(option = 'magma')
plot6 <- ggplot() +
  geom_stars(data = perimeters_r2['DX_class']) +
  #geom_sf(data = perimeters_utm, fill = NA, aes(color = day), lwd = .5) +
  scale_fill_viridis_c(na.value = NA) +
  scale_color_viridis_c(option = 'magma') +
  coord_equal()
plot7 <- ggplot() +
  geom_stars(data = perimeters_r2['day']) +
  geom_raster(data = filter(topocropdf, !is.na(hillshade)), aes(x, y, alpha = hillshade), fill = 'black') +
  scale_fill_viridis_c(option = 'magma', na.value = NA, alpha = .5) +
  scale_alpha_continuous(trans = 'reverse') +
  coord_equal()
ggsave('figures/knp_firedx.pdf', plot5, width = 6, height = 8)
ggsave('figures/knp_firedxclass.pdf', plot6, width = 6, height = 8)
ggsave('figures/knp_firedays.pdf', plot7, width = 6, height = 8)

ggplot() +
  geom_stars(data = perimeters_r2['DX']) +
  geom_raster(data = filter(topocropdf, !is.na(hillshade)), aes(x, y, alpha = hillshade), fill = 'black') +
  scale_fill_viridis_c(na.value = NA, alpha = .8, direction = -1) +
  scale_alpha_continuous(trans = 'reverse') +
  coord_equal() 
ggsave('figures/knp_firedxdem.pdf', width = 6, height = 8)
  
