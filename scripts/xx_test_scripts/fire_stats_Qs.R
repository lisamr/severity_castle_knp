# common stats to put into proposals

library(terra)
library(sf)


veg <- rast('data/SEKI_veg/rSEKI_30m.grd')
lookup <- readRDS('data/SEKI_veg/veg_vars_lookup.rds')
fires <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')
KNP_sev <- rast('outputs/spatial/GEE_CBI/KNP_2021_CBI_bc_corrected.tif')
Cast_sev <- rast('outputs/spatial/GEE_CBI/Castle_2020_CBI_bc_corrected.tif')


names(veg)
names(lookup)
lookup$CalvegDesc %>% sort 


# severity classes?
length(cells(ifel(Cast_sev >= 2.25, 1, NA))) / length(cells(Cast_sev)) # 43.7
cellSize(ifel(Cast_sev >= 2.25, 1, NA), unit = 'ha') %>% 
  global(., 'sum', na.rm = T) #30661.06

# castle + KNP severity
bothfires <- merge(KNP_sev, Cast_sev)
length(cells(ifel(bothfires >= 2.25, 1, NA))) / length(cells(bothfires)) #0.404
cellSize(ifel(bothfires >= 2.25, 1, NA), unit = 'ha') %>% 
  global(., 'sum', na.rm = T) #43052.22


# Qs about mixed conifer
MC <- 30:32
r <- veg$CalVegID
msk <- ifel(r %in% MC, 1, NA)
r_masked <- mask(r, msk) 

# Percentage of mixed conifer forests in the fires:
mask(msk, fires[c(1,3),]) %>% plot(col = 'red')
MC_burned <- mask(msk, fires[c(1,3),]) 
length(cells(MC_burned)) / length(cells(msk)) # 0.427

# of that burned, at high severity?
KNP_Cast <- merge(KNP_sev, Cast_sev) %>% extend(msk)
hi_sev <- ifel(KNP_Cast >= 2.25, 1, NA)
msk_ex <- extend(msk, KNP_Cast)
MC_hisev <- mask(hi_sev, msk_ex) 

length(cells(MC_hisev)) / length(cells(MC_burned)) # 0.362
