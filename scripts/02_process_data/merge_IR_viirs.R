# create sample locations. KNP only for now.
# sample from spatiotemporal data (fire dates, climate)
# use IR data, mask out missing days

rm(list = ls())
library(sf) 
library(terra)
library(tidyverse)
library(tictoc)
library(foreach)
library(lubridate)

# IR + viirs data
IR <- st_read('outputs/spatial/IR/knp_byDay_u.geojson')
viirs_list <- read_rds('outputs/spatial/VIIRS/perims_viirs_suomi.rds') # knp for now




# manipulate IR data a bit ------------------------------------------------

IR <- IR %>% # adjust dates too. flights taken around midnight, dates recorded for following morning
  mutate(date = as_date(date), 
         date_adj = date - 1)


# merge IR and viirs ------------------------------------------------------


# compare growth
p_IR <- IR %>% 
  #filter(date_adj <= '2021-09-12') %>% 
  filter(date_adj <= max(viirs_list$polys_growth$date_adj)) %>% 
  arrange(rev(date_adj)) %>% 
  ggplot() +
  geom_sf( aes(fill = date_adj), color = NA) +
  harrypotter::scale_fill_hp(option = 'Ravenclaw')

p_viirs <- viirs_list$polys_growth %>% 
  #filter(date_adj <= '2021-09-12') %>% 
  ggplot() +
  geom_sf( aes(fill = date_adj), color = NA) +
  harrypotter::scale_fill_hp(option = 'Ravenclaw')

# wow, they are very different...
cowplot::plot_grid(p_IR, p_viirs)

# find missing days for IR
start_date <- as_date('2021-09-11') 
end_date <- as_date('2021-10-08')
date_range <- seq(start_date, end_date, by = 'day')
missing_days <- date_range[!date_range %in% IR$date_adj]

# bind IR with viirs. colnames need to be the same
viirs_missing <- viirs_list$polys_total %>% 
  filter(date_adj %in% missing_days) %>% 
  rename(ha = ha_total, geometry = polygons) %>% 
  mutate(ha = as.numeric(ha), type = 'VIIRS')
perimeters_total <- bind_rows(
  IR %>% select(-c(date, jday)) %>% mutate(type = 'IR'), 
  viirs_missing) %>% 
  arrange(date_adj) %>% 
  mutate(jday_adj = yday(date_adj)) %>% 
  filter(date_adj <= end_date)
print(perimeters_total, n = Inf)

# get growth for each day. i don't want overlapping polygons
f_growth <- function(polys, i){
  suppressWarnings(
    st_difference(polys[i,], st_geometry(polys[i-1,]))
   )  %>% 
     mutate(ha_growth = st_area(.) %>% units::set_units('ha'), .after = ha)
}


# get growths for each day
perimeters_growth <- map_df(2:nrow(perimeters_total), ~ f_growth(perimeters_total, .))

# viirs perimeters overestimate. get intersection between viirs day and next IR day. 
# do it manually
viirs_growth <- perimeters_growth %>% 
  filter(type == 'VIIRS')
perimeters_total %>% print(n = Inf)
growth_0912 <- st_intersection(viirs_growth[1,], perimeters_total[3,]) %>% 
  st_geometry()
growth_1004 <- st_intersection(viirs_growth[2,], perimeters_total[25,]) %>% 
  st_geometry() %>% 
  st_collection_extract('POLYGON') %>% st_union()
growth_1007 <- st_intersection(viirs_growth[3,], perimeters_total[28,]) %>% 
  st_geometry() %>% 
  st_collection_extract('POLYGON') %>% st_union()

# check. there should be differences.
st_area(growth_0912)
st_area(viirs_growth$geometry[1])
st_area(growth_1004)
st_area(viirs_growth$geometry[2])
st_area(growth_1007)
st_area(viirs_growth$geometry[3])


# replace the geometries
perimeters_growth$date_adj
perimeters_growth$geometry[1] <- growth_0912
perimeters_growth$geometry[23] <- growth_1004
perimeters_growth$geometry[26] <- growth_1007

# update the sizes
perimeters_growth <- perimeters_growth %>% 
  mutate(ha_growth = st_area(.) %>% units::set_units('ha'))

perimeters_growth['type'] %>% plot



# export ------------------------------------------------------------------

write_rds(perimeters_growth, 'outputs/spatial/IR/knp_growth_IR_VIIRS.rds')
