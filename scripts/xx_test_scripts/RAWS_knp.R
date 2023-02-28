# figuring out how to work with RAWS data to get info on inversion layers

library(MazamaSpatialUtils)
library(RAWSmet)
library(tidyverse)
library(lubridate)

setRawsDataDir('~/Data/RAWS')
setSpatialDataDir('~/Data/Spatial')
MY_PASSWORD <- 'Fire2010'

# get CA RAWS
#ca_meta <- wrcc_loadMeta(stateCode = 'CA')
#saveRDS(ca_meta, 'outputs/ca_meta_RAWS.rds') #not sure if I need to save, but doing anyways
ca_meta <- read_rds('outputs/ca_meta_RAWS.rds')


# filter down to SEKI only: 
# ash station, cedar grove, mineral king (wolverton), grant grove (park ridge)
# Shadequarter, Case Mountain, Milk Ranch, Sugarloaf, Pinehurst 
unique(ca_meta$siteName) %>% sort
SEKI_stations <- c('Ash Mountain', 'Wolverton', 'Park Ridge', 'Cedar Grove',
                   'Shadequarter', 'Case Mountain', 'Milk Ranch', 'Sugarloaf', 'Pinehurst' )
SEKI_meta <- ca_meta %>% 
  filter(siteName %in% SEKI_stations)

SEKIRAWS_2021 <- wrcc_loadMultiple(
  wrccID = SEKI_meta$wrccID,
  year = 2021, 
  meta = ca_meta,
  password = MY_PASSWORD)

raws_leaflet(SEKI_meta)
raws_leaflet(ca_meta)


# =============================================================================
# OVERVIEW OF DATA
# =============================================================================

# extract the data
SEKIRAWS_2021d <- lapply(SEKIRAWS_2021, function(x) bind_cols(x$meta, x$data)) %>% 
  bind_rows() %>% 
  mutate(datetime = with_tz(datetime, "America/Los_Angeles"), 
         dateround = round_date(datetime, unit = 'day'))


# visualize
names(SEKIRAWS_2021d)
p1 = SEKIRAWS_2021d %>% 
  filter(date(datetime) > '2021-09-01' & date(datetime) < '2021-10-24', 
               hour(datetime) %in% 6:11) %>% 
  group_by( dateround, siteName, elevation) %>% 
  summarise(solarRadiation = mean(solarRadiation)) %>% 
  arrange(dateround) %>% 
  ggplot(., aes(dateround, solarRadiation, group = siteName, color = elevation)) +
  geom_line(lwd = 1) +
  geom_point() +
  scale_color_viridis_c()
p1

p2 = SEKIRAWS_2021d %>% 
  filter(date(datetime) > '2021-09-01' & date(datetime) < '2021-10-15', 
         hour(datetime) %in% 6:11) %>% 
  group_by( dateround, siteName, elevation) %>% 
  summarise(humidity = mean(humidity)) %>% 
  arrange(dateround) %>% 
  ggplot(., aes(dateround, humidity, group = siteName, color = elevation)) +
  geom_line(lwd = 1) +
  geom_point() +
  scale_color_viridis_c(  )
p2

cowplot::plot_grid(p1, p2, nrow = 2)

SEKIRAWS_2021d %>% 
  filter(date(datetime) > '2021-10-01' & date(datetime) < '2021-10-31') %>% 
  group_by( dateround, siteName, elevation) %>% 
  summarise(precipitation = max(precipitation)) %>% 
  arrange(dateround) %>% 
  ggplot(., aes(dateround, precipitation, group = siteName, color = elevation)) +
  geom_line(lwd = 1) +
  geom_point() +
  scale_color_viridis_c(  )




SEKIRAWS_2021d %>% 
  filter(date(datetime) > '2021-09-01' & date(datetime) < '2021-10-15', 
         hour(datetime) %in% 6:11) %>% 
  group_by( dateround, siteName, elevation) %>% 
  summarise(humidity = mean(humidity)) %>% 
  group_by(dateround) %>% 
  mutate(mean = mean(humidity), 
         dHumidity = humidity - mean) %>% 
  ggplot(., aes(dateround, dHumidity, group = siteName, color = elevation)) +
  geom_line(lwd = 1) +
  geom_point() +
  scale_color_viridis_c(  )
