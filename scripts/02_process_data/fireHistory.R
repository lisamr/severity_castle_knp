# prep fire perimeters so that you can get fire history
#(# of fires since 1985, years since last fire up to 1985, severity of last fire)
# severity of last fire is processed in GEE. I read in the outputed shapefiles,
# GEE outputs are saved in the outputs folder. 
library(sf)
library(tidyverse)

# import flattened fire summaries and 2020 and 2021 fire perimeters
frap <- st_read('data/fires/firesummary.shp')
fires <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')
frap <- st_transform(frap, st_crs(fires))


# union frap polygons by a column of interest within each of the fire perimeters. 
f_union <- function(fire, column){
  tmp <- st_intersection(frap, fire)
  tmp %>% 
    # linestrings don't play nice. remove them, keepng only polygons
    st_collection_extract('POLYGON') %>% 
    filter(lsF2019 >= 1985) %>% 
    select(var = {{column}}) %>% 
    group_by(var) %>% 
    summarise(geometry = st_union(geometry)) %>% 
    ungroup 
}

# last fire. send to GEE for getting fire severity for each year's fire.
lastFire_castle <- f_union(fires[1,], lsF2019) %>% rename(lastFire = var)
lastFire_rattlesnake <- f_union(fires[2,], lsF2019) %>% rename(lastFire = var)
lastFire_knp <- f_union(fires[3,], lsF2020) %>% rename(lastFire = var)

# useful to just merge them all together
lastFire_merged <- bind_rows(lastFire_castle, lastFire_rattlesnake, lastFire_knp) %>% 
  group_by(lastFire) %>% 
  summarise(geometry = st_union(geometry)) %>% 
  ungroup %>% 
  mutate(Fire_ID = paste0('Fires_', lastFire),
         Fire_Year = as.integer(lastFire),
         lastFire = NULL)
  

st_write(lastFire_castle, 'outputs/spatial/data4GEE/lastFire_castle.shp', delete_layer = T)
st_write(lastFire_rattlesnake, 'outputs/spatial/data4GEE/lastFire_rattlesnake.shp', delete_layer = T)
st_write(lastFire_knp, 'outputs/spatial/data4GEE/lastFire_knp.shp', delete_layer = T)
st_write(lastFire_merged, 'outputs/spatial/data4GEE/lastFire_merged.shp', delete_layer = T)



# number of fires since 1985.
nFires_castle <- f_union(fires[1,], nF1985_201) %>% rename(nFires85 = var)
nFires_rattlesnake <- f_union(fires[2,], nF1985_201) %>% rename(nFires85 = var)
nFires_knp <- f_union(fires[3,], nF1985_202) %>% rename(nFires85 = var)
plot(nFires_castle)
plot(nFires_rattlesnake)
plot(nFires_knp)
