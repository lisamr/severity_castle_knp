# download weather data from gridmet

# install packages
# remotes::install_github("mikejohnson51/climateR")

# load packages
library(climateR)
library(tidyverse)
library(terra)
library(sf)

# load fires to get aoi
perims <- st_read('data/fires/FRAP_SQF_KNP_perimeters.shp')
aoi_knp <- st_as_sfc(st_bbox(perims[3,]))
aoi_castle <- st_as_sfc(st_bbox(perims[1,]))

# variables for downloading
vars <- c('vpd', 'rmax', 'fm1000', 'tmmx')
dates_knp <- c('2021-09-12', '2021-10-31')
dates_castle <- c('2020-08-23', '2020-11-10')

res_knp <- getGridMET(aoi_knp, varname = vars, 
                   startDate = dates_knp[1], endDate = dates_knp[2])
res_castle <- getGridMET(aoi_castle, varname = vars, 
                      startDate = dates_castle[1], endDate = dates_castle[2])


# export ------------------------------------------------------------------


path <- 'outputs/spatial/weather'
if(!dir.exists(path)) dir.create(path)

walk2(res_knp, vars, 
     ~ writeRaster(.x, glue::glue("{path}/KNP_{.y}.tif"), overwrite = T))
walk2(res_castle, vars, 
      ~ writeRaster(.x, glue::glue("{path}/castle_{.y}.tif"), overwrite = T))



# vis ---------------------------------------------------------------------



# peak at how much temporal variation exists 
quant <- c(.1, .5, .9)
quantile_of_raster <- function(r){
  sapply(r, function(x) global(x, function(y) quantile(y, quant)))
}
tidy_quant <- function(quant_output){
  matrix(unlist(quant_output), nrow = 3) %>% 
    t %>% 
    as_tibble() %>% 
    set_names(paste0('q_', quant)) %>% 
    mutate(day = 1:nrow(.))
} 
qs <- map(res_knp, quantile_of_raster) %>% map(tidy_quant)
qs_knp <- map2_df(qs, names(qs), ~ mutate(.x, var = .y))
qs <- map(res_castle, quantile_of_raster) %>% map(tidy_quant)
qs_castle <- map2_df(qs, names(qs), ~ mutate(.x, var = .y))

bind_rows(qs_knp %>% mutate(fire = 'knp'),
          qs_castle %>% mutate(fire = 'castle')) %>%
  ggplot(., aes(day, q_0.5, color = fire)) +
  geom_pointrange(aes(ymin = q_0.1, ymax = q_0.9), alpha = .5, cex = .4) +
  facet_wrap(~var, scales = 'free_y')
