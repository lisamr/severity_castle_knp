library(tidyverse)
library(sf)
library(lubridate)
library(dbscan)
library(concaveman)
theme_set(theme_classic())

# test out how you'd delineate the perimeters from VIIRS points

# load IR data for comparison
IR <- st_read('outputs/spatial/IR/knp_byDay_u.geojson')

# load VIIRS data
viirs <- st_read('data/VIIRS/SUOMI_VIIRS_C2/fire_archive_SV-C2_352008.shp') %>% 
  st_transform(st_crs(IR))



# get datetime in PDT. time/date currently in UTC.
viirs_d <- viirs %>% 
  mutate(.before = ACQ_DATE, 
         hour = str_sub(ACQ_TIME, 1,2),
         minute = str_sub(ACQ_TIME, 3,4),
         time = paste(hour, minute, sep = ':'),
         datetime = ymd_hm(paste(ACQ_DATE, time)) %>% with_tz('US/Pacific'),
         date = date(datetime),
         jday = yday(date)) %>% 
  select(-c(hour, minute, time, ACQ_DATE, ACQ_TIME)) 


#viirs_d$datetime %>% unique %>% sort

# visualize progress ------------------------------------------------------


distinct(viirs_d, date, jday)
lim <- 260
ggplot() +
  geom_sf(data = IR %>% filter(jday <= lim) %>% arrange(-jday), 
          aes(fill = jday)) +
  geom_sf(data = viirs_d %>% filter(jday <= lim),
          aes(fill = jday), color = 'black', shape = 21) 





# dbscan + concave hull ---------------------------------------------------



f_cluster <- function(pts, d, concavity=2){
  
  # Perform DBSCAN clustering
  dbscan_result <- dbscan(st_coordinates(pts), eps = d, minPts = 3)
  
  # add cluster labels to the original point feature class
  pts$cluster <- dbscan_result$cluster
  
  perims <- pts %>%
    filter(cluster != 0) %>%
    group_split(cluster) %>%
    map(~ concaveman::concaveman(.x, concavity)) %>%
    bind_rows()
  
  return(perims)
}

testdata <- viirs_d 

jdays <- unique(testdata$jday)

polys_list <- testdata %>%
  group_split(jday) %>%
  map(f_cluster, 1125) 

# remove days with no clusters, bind
lengths <- map_vec(polys_list, nrow)
polys <- map2(polys_list[lengths != 0], jdays[lengths != 0], 
     ~ mutate(.x, jday = .y)) %>% 
  bind_rows()





# union polygons from previous day ----------------------------------------



# fire growth should be cumulative. fire perimeter shouldn't shrink.
# you should force polygons to be nested. union with previous day.

f_union <- function(polygons){
  
  # group by day
  polys_byDay <- polygons %>% 
    group_by(jday) %>% 
    summarize()
  
  polys_byDay_u <- polys_byDay
  for (i in 1:(nrow(polys_byDay_u)-1)) {
    poly1 <- polys_byDay_u[i,]
    poly2 <- polys_byDay_u[i + 1,]
    new_geom2 <- st_union(st_geometry(poly1), poly2) 
    
    if(!st_geometry_type(new_geom2) %in% c('MULTIPOLYGON', 'POLYGON')){
      # could make lots of linestrings. keep only polygons
      new_geom2 <- st_collection_extract(new_geom2, 'POLYGON') %>% 
        st_union()
    }
    st_geometry(polys_byDay_u)[i+1] <- new_geom2
    cat('finished unioning jday:', polys_byDay_u$jday[i+1], '\n')
  }
  
  # update area
  polys_byDay_u <- polys_byDay_u %>% 
    mutate(ha_total = st_area(.) %>% units::set_units('ha'), .after = jday)
  
  return(polys_byDay_u)
}

polys_byDay_u <- f_union(polys)

plot_min <- 0
plot_max <- 280
p <- ggplot() +
  geom_sf(data = IR %>% filter(jday <= plot_max, jday >= plot_min) %>% arrange(-jday), 
          aes(fill = jday), alpha = .4) +
  #geom_sf(data = testdata, aes(color = jday)) +
  harrypotter::scale_color_hp() +
  harrypotter::scale_fill_hp() +
  theme(legend.position = 'none')

p + geom_sf(data = polys_byDay_u %>% filter(jday <= plot_max, jday >= plot_min), aes(color = jday), linewidth = 1, fill = NA)
cowplot::plot_grid(
  p + geom_sf(data = polys %>% filter(jday <= plot_max, jday >= plot_min), aes(color = jday), linewidth = 1, fill = NA) ,
    nrow = 1
  
)



# get growth each day -----------------------------------------------------


# find difference so you can see growth each day. polygons shouldn't overlap.

f_growth <- function(polys, i){
  suppressWarnings(
    st_difference(polys[polys$jday == i,], st_geometry(polys[polys$jday == i-1,])) 
    )  %>% 
    mutate(ha_growth = st_area(.) %>% units::set_units('ha'), .after = ha_total)
}
polys_growth <- map_df(polys_byDay_u$jday, ~ f_growth(polys_byDay_u, .x))

p + geom_sf(data = polys_growth %>% filter(jday <= plot_max, jday >= plot_min), 
            aes(color = jday), linewidth = 1, fill = NA) 



full_join(tibble(type_ir = 1, jday = IR$jday), 
          tibble(type_viirs = 1, jday = polys_byDay_u$jday)) %>% 
  left_join(distinct(viirs_d, date, jday)) %>% 
  print(n=Inf)



