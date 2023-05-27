# determine fire perimeters from VIIRS data

library(tidyverse)
library(sf)
library(lubridate)
library(dbscan)
library(concaveman)
theme_set(theme_classic())


# functions ---------------------------------------------------------------

f_cluster_hull <- function(pts, d, concavity=2){
  
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


f_union <- function(polygons, by_daytime){
  
  if(by_daytime){
    # group by day and time of day
    polys_byDay <- polygons %>% 
      group_by(jday_adj, DAYNIGHT) %>% 
      summarize()
  }else{
    # group by day and time of day
    polys_byDay <- polygons %>% 
      group_by(jday_adj) %>% 
      summarize()
  }

  
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
    cat('finished unioning jday_adj:', polys_byDay_u$jday_adj[i+1], '\n')
  }
  
  # update area
  polys_byDay_u <- polys_byDay_u %>% 
    ungroup() %>% 
    mutate(ha_total = st_area(.) %>% units::set_units('ha'), .after = jday_adj )
  
  return(polys_byDay_u)
}

f_growth <- function(polys, i){
  suppressWarnings(
    st_difference(polys[i,], st_geometry(polys[i-1,])) 
  )  %>% 
    mutate(ha_growth = st_area(.) %>% units::set_units('ha'), .after = ha_total)
}

# load data ---------------------------------------------------------------



# load VIIRS data
suomi <- st_read('data/VIIRS/SUOMI_VIIRS_C2/fire_archive_SV-C2_352008.shp') %>% 
  st_transform(crs = 'epsg:32611')
j1 <- st_read('data/VIIRS/J1_VIIRS_C1/fire_nrt_J1V-C2_352007.shp') %>% 
  st_transform(crs = 'epsg:32611')
viirs <- bind_rows(suomi, j1)

# get datetime in PDT. time/date currently in UTC.
viirs_d <- viirs %>% 
  mutate(.before = ACQ_DATE, 
         hour = str_sub(ACQ_TIME, 1,2),
         minute = str_sub(ACQ_TIME, 3,4),
         time = paste(hour, minute, sep = ':'),
         datetime = ymd_hm(paste(ACQ_DATE, time)) %>% with_tz('US/Pacific'),
         date = date(datetime),
         date_adj = case_when(DAYNIGHT == 'N'~date - 1, T~date),
         jday = yday(date),
         jday_adj = ifelse(DAYNIGHT == 'N', jday - 1, jday),
         jday_daynight_adj = paste(jday_adj, DAYNIGHT, sep = '_')) %>% 
  select(-c(hour, minute, time, ACQ_DATE, ACQ_TIME)) %>% 
  filter(CONFIDENCE != 'l') # remove low confidence


viirs_d %>% 
  mutate(hr = hour(datetime)) %>% 
  count(DAYNIGHT, hr)
# D = 12-15pm, N = 1-3am

# make perimeters ---------------------------------------------------------



# distance threshold
d <- 1125

make_perimeters <- function(pts, by_daytime){
  
  # ======================
  # dbscan + concave hull
  # ======================
  
  message('startying to make concave hulls')
  
  if(by_daytime){
    # make a list of polygon hulls
    polys_list <- pts %>%
      group_split(jday_adj, DAYNIGHT) %>%
      map(f_cluster_hull, d) 
    groups <- pts %>%
      group_split(jday_adj, DAYNIGHT) %>% 
      map_df(~ distinct(.x, jday_adj, DAYNIGHT))
    # add in dates
    dates <- pts %>% 
      distinct(jday_adj, date_adj, DAYNIGHT)
  }else{
    # make a list of polygon hulls
    polys_list <- pts %>%
      group_split(jday_adj) %>%
      map(f_cluster_hull, d) 
    groups <- pts %>%
      group_split(jday_adj) %>% 
      map_df(~ distinct(.x, jday_adj))
    # add in dates
    dates <- pts %>% 
      distinct(jday_adj, date_adj)
  }
  
  
  # remove days with no clusters, bind
  lengths <- map_vec(polys_list, nrow)
  if(length(lengths) != nrow(groups)){
    stop( 'polygon list isnt same length as your grouping variables')
  }
  groups <- groups[lengths != 0,]
  polys <- map2(polys_list[lengths != 0], 1:nrow(groups), 
                ~ mutate(.x, groups[.y,])) %>% 
    bind_rows()
  
  # ======================
  # union polygons from previous time step 
  # ======================
  
  message('unioning polygons from previous time step')
  # fire growth should be cumulative. fire perimeter shouldn't shrink.
  # you should force polygons to be nested. union with previous day.
  polys_bygroup_u <- f_union(polys, by_daytime)
  
  # ======================
  # get growth each day 
  # ======================
  
  message('getting growth of each time step')
  # find difference so you can see growth each day. polygons shouldn't overlap.
  polys_growth <- bind_rows(
    polys_bygroup_u[1,] %>% mutate(ha_growth = NA, .after = ha_total),
    map_df(1:nrow(polys_bygroup_u), 
           ~ f_growth(polys_bygroup_u, .x))
  ) 

  
  return(list(polys_total = polys_bygroup_u %>% left_join(dates),
              polys_growth = polys_growth %>% left_join(dates)))
}

perims_viirs <- make_perimeters(viirs_d, F)
perims_viirs_bydaytime <- make_perimeters(viirs_d, T)
perims_viirs_suomi <- make_perimeters(viirs_d %>% filter(SATELLITE == 'N'), F)
perims_viirs_suomi_bydaytime <- make_perimeters(viirs_d %>% filter(SATELLITE == 'N'), T)




# export ------------------------------------------------------------------

path <- 'outputs/spatial/VIIRS/'
if(!dir.exists(path)) dir.create(path)
write_rds(perims_viirs, glue::glue({path},"perims_viirs.rds"))
write_rds(perims_viirs_bydaytime, glue::glue({path},"perims_viirs_bydaytime.rds"))
write_rds(perims_viirs_suomi, glue::glue({path},"perims_viirs_suomi.rds"))
write_rds(perims_viirs_suomi_bydaytime, glue::glue({path},"perims_viirs_suomi_bydaytime.rds"))



# visualize ---------------------------------------------------------------
# 
# perims_viirs$polys_growth
# perims_viirs_bydaytime$polys_growth
# cowplot::plot_grid(
#   ggplot(perims_viirs_bydaytime$polys_growth, aes(fill = DAYNIGHT)) +
#     geom_sf(),
#   ggplot(perims_viirs_suomi_bydaytime$polys_growth, aes(fill = DAYNIGHT)) +
#     geom_sf()
# )

# IR <- st_read('outputs/spatial/IR/knp_byDay_u.geojson')
# IR_growth <- map_df(2:nrow(IR), 
#        ~ suppressWarnings(
#            st_difference(IR[.x,], st_geometry(IR[.x-1,])) 
#          ))
# 
# plot_min <- 0
# plot_max <- 280
# p <- ggplot() +
#   geom_sf(data = IR_growth %>% filter(jday <= plot_max, jday >= plot_min) %>% arrange(-jday), 
#           aes(fill = jday), alpha = .4) +
#   harrypotter::scale_color_hp() +
#   harrypotter::scale_fill_hp() +
#   theme(legend.position = 'none')
# 
# 
# pgrowth <- p + geom_sf(data = perims_viirs$polys_total %>% 
#                                filter(jday <= plot_max, jday >= plot_min),
#                             aes(color = jday, fill = jday), linewidth = 1, fill = NA)
# pgrowth_bydaytime <- p + geom_sf(data = perims_viirs_bydaytime$polys_total %>% 
#                               filter(jday <= plot_max, jday >= plot_min),
#                             aes(color = jday), linewidth = 1, fill = NA)
# 
# 
# cowplot::plot_grid(
#   pgrowth,
#   pgrowth_bydaytime
# )




