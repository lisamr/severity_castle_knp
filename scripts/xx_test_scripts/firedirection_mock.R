# how can you declare if a fire is heading, backing or flanking?
rm(list=ls())
library(sf)
library(tidyverse)
library(fasterize)
library(nngeo) #remotes::install_github("michaeldorman/nngeo")
library(cowplot)
library(raster)
theme_set(theme_classic())



#=========================================================================
# CREATE POLYGONS
#=========================================================================

# load topogrpahy
topo <- brick('../../../Documents/ArcGIS/Spatial_data/from_GEE/topography30m.tif')
randompt <- as.data.frame(sampleRandom(topo, 1, xy = T))
X0 = randompt$x
Y0 = randompt$y

# make polygons, where one is inside the other
p3 <- st_buffer(st_point(c(X0 + 800,Y0)), 1000)
p2 <- st_buffer(st_point(c(X0 + 400,Y0)), 500)
p1 <- st_buffer(st_point(c(X0,Y0)), 100)
geometry <- st_sfc(p1, p2, p3)
ndays <- length(geometry)
polygons <- st_sf(day = 1:ndays, geometry, crs = crs(topo)) 

# crop topo
topocrop <- crop(x = topo, polygons)
plot(topocrop, col = grey.colors(20))

#=========================================================================
# TURN POLYGONS INTO RASTER
#=========================================================================

# convert to raster
pr <- fasterize(polygons, topocrop$elevation, field = 'day', fun = 'min')
names(pr) <- 'day'

plot(topocrop$hillshade, col = grey.colors(10))
#plot(pr)
plot(polygons$geometry, add = T)



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
pr_pts <- rasterToPoints(pr, spatial = T) %>% st_as_sf() %>% arrange(day) #raster as points
perim_pts <- st_cast(polygons, 'LINESTRING') %>%  #perimeters as lines, ignore warning?
  st_line_sample(., density= 5) %>% #convert to multipoints
  st_sf(day = 1:ndays, .) %>% #add attributes
  st_cast('POINT') #convert to points
perim_pts$ID <- 1:nrow(perim_pts) #add ID to perimeter points for keeping track

# visualize whats going on
plot(pr)
plot(pr_pts, add = T) 
plot(perim_pts[,'day'], add = T, col = 'black') 


#find the nearest edge and add to sf object
nearestedge <- rep(NA, nrow(pr_pts))
v1 <- 1
for(i in 1:ndays){
  index <- st_nearest_feature(filter(pr_pts, day == i), filter(perim_pts, day == i-1))
  ID <- filter(perim_pts, day == i-1)[index,] %>% pull('ID')
  v2 <- length(index) -1 + v1
  nearestedge[v1:v2] <- ID 
  v1 <- v2 + 1
}
pr_pts$nearestedge <- nearestedge
nearestedge1 <- nearestedge
nearestedge1[is.na(nearestedge1)] <- 1 #st_azimuth doesn't allow NAs, so temporarily assign NA values as 1
pr_pts$nearestedge1 <- nearestedge1


#calculate azimuth (and distance so I know I'm doing it right) to nearest edge
pr_pts$dist <- st_distance(pr_pts, perim_pts[nearestedge,], by_element = T) %>% as.numeric()
pr_pts$az <- st_azimuth(pr_pts, perim_pts[nearestedge1,])
pr_pts$az[is.na(nearestedge)] <- NA #change the values within first polygon as NA


# viz. I think that looks right.
plot1 <- ggplot() +
  geom_sf(data = pr_pts, aes(color = az), size = 3) +
  scale_color_viridis_c(na.value = NA)
plot2 <- ggplot() +
  geom_sf(data = pr_pts, aes(color = dist), size = 3) +
  scale_color_viridis_c(na.value = NA)
plot3 <- ggplot() +
  geom_sf(data = pr_pts, aes(color = nearestedge), size = 3) +
  scale_color_viridis_c(na.value = NA)
plot4 <- ggplot() +
  geom_sf(data = pr_pts, aes(color = day), size = 3) +
  geom_sf(data = perim_pts, size = 1) +
  scale_color_viridis_c(na.value = NA)
plot_grid(plot1, plot2, plot3, plot4)

ggplot() +
  geom_sf(data = perim_pts, aes(color = ID)) +
  scale_color_viridis_c(na.value = NA)




#=========================================================================
# DIRECTION OF FIRE
#=========================================================================

# first apply smoothing filter on aspect.
topocrop$aspect_sm <- focal(topocrop$aspect, w=matrix(1/9,nrow=3,ncol=3)) 

# get topo values for all the raster points
topovals <- extract(topocrop, pr_pts)
pr2 <- cbind(pr_pts, topovals)

# get fire direction
pr2 <- pr2 %>%
  mutate(DX = abs(az - aspect),
         DX = ifelse(DX>180, 360-DX, DX),
         DX_class = case_when(
           DX <= 60 ~ 'H',
           DX <= 120 ~ 'F',
           T ~ 'B'
         ),
         DX_class = ifelse(is.na(DX), NA, DX_class))
plot5 <- ggplot() +
  geom_sf(data = pr2, aes(color = DX), shape = 15) +
  geom_sf(data = polygons, fill = NA) +
  scale_color_viridis_c(direction = -1)
plot5b <- ggplot() +
  geom_sf(data = pr2, aes(color = DX_class), shape = 15) +
  geom_sf(data = polygons, fill = NA) +
  scale_color_viridis_d()
plot6 <- ggplot() +
  geom_sf(data = pr2, aes(color = hillshade), shape = 15) +
  geom_sf(data = polygons, fill = NA) +
  scale_color_gradient(low = grey(.1), high = grey(.9))
#plot6
plot_grid(plot5, plot5b, plot6, nrow = 1)


#=========================================================================
# TRASH
#=========================================================================


# create a raster layer of the perimeter of each day, 
# create a distance matrix between all the pixels of one day and the perimeter,
# identify the nearest perimeter pixel,
# calculate the azimuth to it.

# get perimeters, keep as a Raster object to calculate distances
plines <- st_cast(polygons, 'MULTILINESTRING')
prlines <- stars::st_rasterize(plines, stars::st_as_stars(r), options = 'All_TOUCHED=T') %>% 
  as('Raster')
plot(prlines)
plot(plines, add = T)





# get distance matrix between all pixels and pixels from day before
prlines
pr

head(pr)
head(prlines)
plot(pr)
plot(prlines, add = T, col = c('black', 'red'))
xy2 <- raster::xyFromCell(pr, which(pr[] == 2)) # coordinates of day two
xy1_line <- raster::xyFromCell(prlines, which(prlines[] == 1)) # coordinates of day 1 boundary

dist()
dist(xy2[1,], xy1_line)


