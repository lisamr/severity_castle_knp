library(MazamaSpatialUtils)
library(PWFSLSmoke)
dir.create('~/Data/Spatial', recursive=TRUE)
setSpatialDataDir('~/Data/Spatial')
installSpatialData()

MKrd <- 'lon_.118.611_lat_36.453_apcd.1032'
AshMt <- '061070009_01'
TR <- 'lon_.118.912_lat_36.428_apcd.1021'
HillTop <- 'lon_.118.962_lat_36.739_apcd.1033'


knp_fire <-
  monitor_loadAnnual(2021) %>%
  monitor_subset(stateCodes = 'CA') %>%
  monitor_subset(tlim = c(20210901, 20211015))

monitor_leaflet(knp_fire)


monitor_timeseriesPlot(
  Sacramento,
  style='aqidots',
  pch=16,
  xlab="2018"
)

knp_sites <- knp_fire %>% 
  monitor_subset(monitorIDs = c(TR))
monitor_timeseriesPlot(
  knp_sites, style='aqidots', pch = 16
)
