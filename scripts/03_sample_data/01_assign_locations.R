# sample points from within fire perimeter. Constraints:
# - lie within area where 90% of the growth occurred. determine from IR. 
# - on a grid? distance between points >= 2*ring_radius_max


library(tidyverse)
library(sf)

# load things -------------------------------------------------------------


IR <- st_read('outputs/spatial/IR/knp_byDay_u.geojson')



# =================
# set arguments
# =================

N <- 100 # N samples
growth_threshold <- .99 # threshold for cumulative growth
max_radius <- 1500 # maximum radius of the ring



# define sampling area ----------------------------------------------------

# dont want to sample in areas where growth is minimal. 
perc_ha <- IR %>% 
  mutate(p_ha = ha/max(ha)) %>% 
  pull(p_ha)
last_day <- max(IR$date[perc_ha < growth_threshold])
fire_perimeter <- IR %>% 
  filter(date == last_day) 



# define sampling grid ----------------------------------------------------

# crap, i'm gonna need to deal with spatial autocorrelation. 
# a grid with 2000m between samples is not gonna cut it. ~150 samples?!?

tmp <- st_make_grid(fire_perimeter, square = F, cellsize = 2000)
sqrt(st_area(tmp))
length(tmp)

tmp_inter <- st_intersection(tmp, fire_perimeter)

ggplot() +
  geom_sf(data = fire_perimeter) +
  geom_sf(data = tmp_inter, 
          fill = NA, lwd = .5)

as.numeric(st_area(fire_perimeter) / 900) / 4
