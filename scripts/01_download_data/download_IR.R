# IR Data. assume the date is most likely midnight (early morning or evening of previous day)

library(glue)
library(sf)
library(RCurl)
library(XML)
library(tidyverse)
library(fs)
library(lubridate)

source('scripts/functions/functions_download_from_URL.R')
useCRS <- 32611 # need to transform since not all consistent. wgs84/11N

# load fire perimeters from FRAP. IR data mixes SQF and Rattlesnake halfway through.
# need to either seperate or just mask out. 
fires <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')

# download IR data --------------------------------------------------------


# KNP SHAPEFILES ==================================

# Set the URL of the parent directory
url <- "https://ftp.wildfire.gov/public/incident_specific_data/calif_s/2021_Incidents/CA-KNP-000122_KNPComplex/IR/"

# get paths to the gdb files
subfolders <- f_get_subfolders(url)

# Download the file from the URL and save it to the local path
dwn_directory <- 'data/IR/shp/KNP'
if (!dir_exists(dwn_directory)) dir_create(dwn_directory, recurse = T)
shp_urls_knp <- map(subfolders, f_get_file_paths_r, 'hapefiles.zip', c('Colony', 'Paradise'))
# manually add missing links that don't have the suffix
shp_urls_knp[[22]] <- 'https://ftp.wildfire.gov/public/incident_specific_data/calif_s/2021_Incidents/CA-KNP-000122_KNPComplex/IR/20211006/20211006_KNP%20Complex_IR%20Polygon.zip'
shp_urls_knp <- unlist(shp_urls_knp)

# download files, unzip file, deletes zipped file
walk(shp_urls_knp, ~ f_download(.x, dwn_directory, T))

# read all shp files as sf object
shp_paths_knp <- dir_ls('data/IR/shp/KNP', recurse = T, glob = '*HeatPerimeter.shp$|*Polygon.shp')
knp_perims_shp <- map(shp_paths_knp, st_read)
knp_perims_shp <- map(knp_perims_shp, ~ st_transform(.x, crs = useCRS))



# SQF =================================

# repeat for SQF fire
url <- 'https://ftp.wildfire.gov/public/incident_specific_data/calif_s/2020_Incidents/CA-SQF-002622_SQF%20COMPLEX/IR/NIROPS/'

# they uploaded .kmz and shapefiles for this fire. do both and cross reference.
# get urls to the files
subfolders <- f_get_subfolders(url)
shp_urls_sqf <- map(subfolders, ~ f_get_file_paths(.x, suffix = '.zip', exclude = 'Products'))
shp_urls_sqf <- unlist(shp_urls_sqf)

# download the files, unzip into own folder, delete zipped files
dwn_directory_shp <- 'data/IR/shp/SQF'
if (!dir_exists(dwn_directory_shp)) dir_create(dwn_directory_shp, recurse = T)
walk(shp_urls_sqf, ~ f_download(.x, dwn_directory_shp, T))

# read all files as sf object
shp_paths_sqf <- dir_ls('data/IR/shp/SQF', recurse = T, glob = '*Perimeter*.shp$')
# keep only files with perimeter in basename
shp_paths_sqf <- shp_paths_sqf[grep('[Pp]erimeter', basename(shp_paths_sqf))] 
length(shp_paths_sqf)

# read in shp. some files are corrupt with the zplane, use readOGR, then convert to sf.
sqf_perims_sp <- map(shp_paths_sqf, rgdal::readOGR) 
sqf_perims_sf <- map(sqf_perims_sp, st_as_sf) %>% map(~ st_transform(.x, crs = useCRS))




# group polygons by day and aggregate -------------------------------------

# KNP =============================

# not sure which dates to use. file name, createdate, datecurrent, or polygondat?
# after inspecting it, shapefile dates got truncated, but rounded by day. use createdate,
# which i think is closest to the date/time on file name. seems like file name
# was manually entered, so rely on the automaticly generated fields.

# look at the shapefiles. 
knp_all <- knp_perims_shp %>% 
  bind_rows() %>% 
  mutate(fire = 'KNP', 
         date = ymd(CreateDate)) %>% 
  filter(FeatureCat == 'IR Heat Perimeter') %>% 
  select(fire, date)

# summarize by day
knp_byDay <- knp_all %>% 
  st_make_valid() %>% # had issues with duplicated vertices
  group_by(date) %>% 
  summarize() %>% 
  mutate(ha = as.numeric(st_area(.)) * .0001, 
         jday = yday(date),
         .after = date)



# SQF =======================================================

# columns aren't consistent across days.
# not all the files have dates either. Get dates from file name. usually in the 
# ymd_hs format. need to fix a few manually. 
filenames_sqf <- basename(names(sqf_perims_sf))
filenames_sqf[85:90] <- gsub('.shp', '_0000.shp', filenames_sqf[85:90]) # rounds to nearest day
date_strings <- str_extract(filenames_sqf, '\\d{8}_c\\d{4}|\\d{8}_\\d{4}') %>% 
  gsub('c', '', .) %>% 
  ymd_hm()
date_strings[is.na(date_strings)] <- ymd_hm('20200912_2200')

# round to nearest day
sqf_dates_best <- round_date(date_strings, 'day') %>% ymd

# delete 9/6. after digging, realized Castle was missed and it's only the shotgun fire. 
idx_delete <- which(sqf_dates_best == '2020-09-06')
sqf_dates_best <- sqf_dates_best[-idx_delete] 
sqf_perims_sf <- sqf_perims_sf[-idx_delete]

# finally, bind and group by day
sqf_all <- map2(sqf_perims_sf, sqf_dates_best, ~ transmute(.x, date = .y)) %>% 
  bind_rows()
# summarize by day
sqf_byDay <- sqf_all %>% 
  group_by(date) %>% 
  summarize() %>% 
  mutate(ha = as.numeric(st_area(.)) * .0001, 
         jday = yday(date),
         .after = date)





# ensure matches up with frap perimeters ----------------------------------

# # Castle
# ggplot() +
#   geom_sf(data = sqf_byDay[nrow(sqf_byDay),], fill = 'blue', alpha = .3) +
#   geom_sf(data = fires[1,], fill = NA, color = 'red')
# 
# # rattlesnake
# ggplot() +
#   geom_sf(data = sqf_byDay[nrow(sqf_byDay),], fill = 'blue', alpha = .3) +
#   geom_sf(data = fires[2,], fill = NA, color = 'red')

# plan...remove rattlesnake from the analysis. better to choose where it intersects 
# the castle fire. it won't leave a ring around rattlesnake and you also dont 
# have complete information on shotgun fire. 
castle_byDay <- st_intersection(sqf_byDay, fires$geometry[fires$FIRE_NAME == 'CASTLE']) %>% 
  # recalculate area per day
  mutate(ha = as.numeric(st_area(.)) * .0001)

# do the same for KNP fire
# KNP looks fine
knp_byDay <- st_intersection(knp_byDay, 
                             fires$geometry[fires$FIRE_NAME == 'KNP Complex']) %>% 
  mutate(ha = as.numeric(st_area(.)) * .0001) # recalculate area per day


# plot
ggplot() +
  geom_sf(data = castle_byDay[nrow(castle_byDay),], fill = 'blue', alpha = .3) +
  geom_sf(data = fires[1,], fill = NA, color = 'red')

ggplot() +
  geom_sf(data = knp_byDay[nrow(knp_byDay),], fill = 'blue', alpha = .3) +
  geom_sf(data = fires[3,], fill = NA, color = 'red')




# find discrepancies by day -----------------------------------------------


# castle fire ================================================

# can't do geometry calculations on geometry collections. 
geomtypes <- st_geometry_type(castle_byDay)

# pull out just polygons from geometry collections
fix_gc <- castle_byDay[geomtypes == 'GEOMETRYCOLLECTION',]
castle_byDay$geometry[which(geomtypes == 'GEOMETRYCOLLECTION')] <- 
  st_collection_extract(fix_gc, 'POLYGON') %>% 
  st_union()

# each day should encompass the previous day. ok if minor differences, but 
# shouldn't be huge. this would indicate polygon is missing chunks.
res <- rep(NA, nrow(castle_byDay)-1)
for(i in 1:length(res)){
  poly1 <- castle_byDay[i,]
  poly2 <- castle_byDay[i + 1,]
  res[i] <- st_difference(poly1, poly2) %>% st_area
  cat('finished i:', i, '\n')
}
castle_diff <- as_tibble(castle_byDay) %>% 
  mutate(growth = c(0, diff(ha)),
         missing = c(0, round(res*.0001, 1)),
         geometry = NULL)
print(castle_diff, n = Inf)

# you should force polygons to be nested. union with previous day.
castle_byDay_u <- castle_byDay
for (i in 1:(nrow(castle_byDay_u)-1)) {
  poly1 <- castle_byDay_u[i,]
  poly2 <- castle_byDay_u[i + 1,]
  new_geom2 <- st_union(poly1$geometry, poly2) 
  # makes lots of linestrings. keep only polygons
  new_geom2 <- st_collection_extract(new_geom2, 'POLYGON') %>% st_union()
  castle_byDay_u$geometry[i + 1] <- new_geom2
  cat('finished i:', i, '\n')
}
# update area
castle_byDay_u$ha <- as.numeric(st_area(castle_byDay_u)) * .0001

# checks
castle_diff$diff_union <- round(castle_byDay_u$ha - castle_byDay$ha, 2)
print(castle_diff, n = Inf)
diff(castle_byDay_u$ha)


# knp fire ================================================

# can't do geometry calculations on geometry collections. its all good
st_geometry_type(knp_byDay)

# each day should encompass the previous day. ok if minor differences, but 
# shouldn't be huge. this would indicate polygon is missing chunks.
res <- rep(NA, nrow(knp_byDay)-1)
for(i in 1:length(res)){
  poly1 <- knp_byDay[i,]
  poly2 <- knp_byDay[i + 1,]
  res[i] <- st_difference(poly1, poly2) %>% st_area
  cat('finished i:', i, '\n')
}
knp_diff <- as_tibble(knp_byDay) %>% 
  mutate(growth = c(0, diff(ha)),
         missing = c(0, round(res*.0001, 1)),
         geometry = NULL)
print(knp_diff, n = Inf)

# you should force polygons to be nested. union with previous day.
knp_byDay_u <- knp_byDay
for (i in 1:(nrow(knp_byDay_u)-1)) {
  poly1 <- knp_byDay_u[i,]
  poly2 <- knp_byDay_u[i + 1,]
  new_geom2 <- st_union(poly1$geometry, poly2) 
  # makes lots of linestrings. keep only polygons
  new_geom2 <- st_collection_extract(new_geom2, 'POLYGON') %>% st_union()
  knp_byDay_u$geometry[i + 1] <- new_geom2
  cat('finished i:', i, '\n')
}
knp_byDay_u$ha <- as.numeric(st_area(knp_byDay_u)) * .0001


# plot it
knp_byDay_u %>% 
  arrange(rev(jday)) %>% 
  ggplot() +
  geom_sf(aes(fill = jday)) +
  scale_fill_viridis_c()
as_tibble(knp_byDay) %>% 
  transmute(date, jday, ha, diff = c(0, diff(ha))) %>% 
  mutate(ha_u = knp_byDay_u$ha,
         diff_u = c(0, diff(knp_byDay_u$ha))) %>% 
  print(n = Inf)
st_geometry_type(knp_byDay_u)



# export data -------------------------------------------------------------

if(!dir_exists('outputs/IR')) dir_create('outputs/IR')
st_write(castle_byDay, 'outputs/IR/castle_byDay.geojson', delete_dsn = T)
st_write(castle_byDay_u, 'outputs/IR/castle_byDay_u.geojson', delete_dsn = T)
st_write(knp_byDay, 'outputs/IR/knp_byDay.geojson', delete_dsn = T)
st_write(knp_byDay_u, 'outputs/IR/knp_byDay_u.geojson', delete_dsn = T)





# trash -------------------------------------------------------------------


# 
# # # KNP GDB FILES ==================================
# # # Set the URL of the parent directory
# url <- "https://ftp.wildfire.gov/public/incident_specific_data/calif_s/2021_Incidents/CA-KNP-000122_KNPComplex/IR/"
# 
# # get paths to the gdb files
# subfolders <- f_get_subfolders(url)
# gdb_paths <- map(subfolders, f_get_file_paths_r, '.gdb.zip', c('Colony', 'Paradise')) %>%
#   unlist
# 
# # Download the file from the URL and save it to the local path
# dwn_directory <- 'data/IR/gdb/KNP'
# if (!dir_exists(dwn_directory)) dir_create(dwn_directory, recurse = T)
# walk(gdb_paths, ~ f_download(.x, dwn_directory, F))
# 
# # read as sf object
# knp_perims_gdb <- map(file.path(dwn_directory, basename(gdb_paths)),
#     ~ st_read(.x, layer = 'IR_Polygon'))
# 
# # #KNP =====================================================
# knp_all_gdb <- knp_perims_gdb %>%
#   bind_rows() %>%
#   filter(FeatureCategory == 'IR Heat Perimeter') %>% 
#   transmute(fire = 'KNP',
#          polygonDate = ymd_hms(PolygonDateTime),
#          CreateDate = ymd_hms(CreateDate),
#          DateCurrent  = ymd_hms(DateCurrent )
#          ) 
# 
# # time of file corresponds to create date
# names(knp_perims_shp)[3]
# gdb_paths[3]
# knp_perims_gdb[[3]]
# knp_perims_shp[[3]]
# 
# st_geometry_type(knp_all_gdb) =='MULTISURFACE' %>% unique
#  
# knp_all_gdb %>% 
#   filter(st_geometry_type(.) == 'MULTISURFACE')
# 
# 
# knp_all_gdb %>% st_area() # wont work, cant figure it out
# 

# SQF ========================================================
# 
# kmz_urls <- map(subfolders, ~ f_get_file_paths(.x, suffix = '.kmz|.kml'))
# kmz_urls <- unlist(kmz_urls)
# dwn_directory_kmz <- 'data/IR/kmz/SQF'
# if (!dir_exists(dwn_directory_kmz)) dir_create(dwn_directory_kmz, recurse = T)
# walk(kmz_urls, ~ f_download(.x, dwn_directory_kmz, T))
# # needs to be manually added
# download.file('https://ftp.wildfire.gov/public/incident_specific_data/calif_s/2020_Incidents/CA-SQF-002622_SQF%20COMPLEX/IR/NIROPS/20200928/20200928_SQF.kml',
#               file.path(dwn_directory_kmz, '20200928_SQF', '20200928_SQF.kml'))
# kml_paths <- dir_ls('data/IR/kmz/SQF', recurse = T, glob = '*.kml')
# # need to figure out which layer to read...
# noPerim <- map(kml_paths, st_layers) %>%  # choose layer name *containing* "Perimeter"
#   map(~ grep('Perimeter', .x$name, value = F)) 
# 
# # which ones don't have a single match to 'Perimeter'?
# trouble_kmls <- tibble(file = kml_paths, 
#                        perims = map_vec(noPerim, length) %>% unname,
#                        index = seq_along(kml_paths)
# ) %>% 
#   filter(perims != 1)
# # look at them. do a case_when statement. 
# map(kml_paths[trouble_kmls$index], st_layers) %>% 
#   map(~ .x$name)
# f_select_layer <- function(i){
#   lyrs <- st_layers(i)$name
#   perims <- grep('Perimeter', lyrs, value = T)
#   val <- paste0(lyrs[!lyrs == 'Isolated Heat Sources'], collapse = ' ') 
#   pLen <- length(perims)
#   layer <- case_when(pLen == 1 ~ perims[1],
#                      pLen == 2 ~ 'Fire Perimeter',
#                      pLen == 0 & length(lyrs) == 2 ~ val,
#                      T ~ NA
#   )
#   return(layer)
# }
# 
# # extracts first string of 8 digits.
# kml_dates <- str_extract(kml_paths, "\\d{8}") %>%  
#   lubridate::ymd(.)
# 
# # read in files
# kml_layers <- map(kml_paths, f_select_layer)
# sqf_perims <- pmap(.l = list(kml_paths,
#                              kml_layers,
#                              kml_dates),
#                    .f = ~ st_read(..1, layer = ..2) %>% 
#                      select(-Description) %>% 
#                      mutate(date = ..3)
# )
# sqf_perims <- unname(sqf_perims)
# 
# 
# 
# 
# 
# 
# # SQF--clean up, then merge. following need to be filtered to include only the heat perimeters. 
# idx <- trouble_kmls %>% 
#   filter(perims == 0) %>% 
#   pull(index)
# sqf_perims[idx] <- map(sqf_perims[idx], ~ filter(.x, Name == 'Heat Perimeter'))
# # now merge
# sqf_all <- bind_rows(sqf_perims) %>% 
#   mutate(Name = NULL, fire = 'SQF', jDay = yday(date)) %>% 
#   arrange(rev(date))
# dir_create('outputs/IR/SQF')
# st_write(sqf_all, 'outputs/IR/SQF/sqf_IR.geojson', delete_dsn = T)
# 
# 
# tmp <- sqf_all %>% 
#   mutate(area = st_area(.), ha = as.numeric(area) * 0.0001 ) 
# as_tibble(tmp) %>% 
#   group_by(date) %>% 
#   summarise(ha = sum(ha)) %>% 
#   mutate(diff = c(0, diff(ha))) %>% 
#   print(n = Inf)
# 
# # this IR data includes the rattlesnake fire and I think that's where the pulse 
# # on 10-14 came from.
# sqf_all %>% 
#   filter(date != ymd('2020-10-14')) %>% 
#   ggplot(.) +
#   geom_sf(aes(fill = date)) +
#   scale_fill_viridis_c()
# 
# 
# cowplot::plot_grid(
#   sqf_all %>% 
#     filter(date == ymd('2020-10-13')) %>% 
#     #filter(date >= ymd('2020-10-14') & date <= ymd('2020-10-15')) %>% 
#     arrange(rev(date)) %>% 
#     ggplot() +
#     geom_sf() +
#     labs(title = '13'),
#   sqf_all %>% 
#     filter(date == ymd('2020-10-14')) %>% 
#     #filter(date >= ymd('2020-10-14') & date <= ymd('2020-10-15')) %>% 
#     arrange(rev(date)) %>% 
#     ggplot() +
#     geom_sf() +
#     labs(title = '14'),
#   sqf_all %>% 
#     filter(date == ymd('2020-10-15')) %>% 
#     #filter(date >= ymd('2020-10-14') & date <= ymd('2020-10-15')) %>% 
#     arrange(rev(date)) %>% 
#     ggplot() +
#     geom_sf() +
#     labs(title = '15')
# )
