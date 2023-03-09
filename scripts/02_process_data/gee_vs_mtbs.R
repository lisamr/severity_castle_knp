# cbi data from GEE seems really off. ~3% of knp > 2.25. 
# compare to mtbs

library(terra)
library(tidyverse)


# load data
cbi <- rast('outputs/spatial/GEE_CBI/KNP_2021_CBI_bc.tif')
mtbs <- rast('data/fires/mtbs_SQF_KNP.grd')
comp <- rast('../master_spatial/spatial_data/GEE_composite_severity/KNP_COMPLEX_rdnbr.tif')

#this function was used to convert rdnbr to cbi and ba
extendedAssessments <- function(RdNBR){
  RdNBR[RdNBR <= -369] <- -368.99999999
  CBI <- 1/.388 * log((RdNBR + 369)/421.7)
  CBI[CBI < 0] <- 0
  CBI[CBI > 3] <- 3

  # need to deal with extreme values.
  RdNBR[RdNBR < 166.5] <- 166.5
  RdNBR[RdNBR > 389*2] <- 389*2
  BA <- (sin((RdNBR - 166.5)/389))^2 * 100
  br <- rast(list(RdNBR, CBI, BA))
  names(br) <- c('RdNBR', 'CBI', 'BA')
  return(br)
}
comp <- extendedAssessments(comp)

# force rasters to common extent. use gee layer.
mtbs_c <- crop(mtbs, cbi)
comp_c <- crop(comp, cbi) %>% extend(cbi)
brick <- list(cbi, comp_c) %>% rast

# compare
plot(brick[[c('CBI_bc', 'CBI')]])
brick_pts <- spatSample(brick, 5000, na.rm = T)
plot(brick_pts$CBI, brick_pts$CBI_bc, ylab = 'GEE cbi', xlab = 'ref cbi', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

df <- data.frame(gee_cbi = brick_pts$CBI_bc, 
           ref_cbi = brick_pts$CBI)

cowplot::plot_grid(
  ggplot(df, aes(ref_cbi, gee_cbi)) +
    geom_point(alpha = .2) +
    geom_abline(slope = 1, lwd = 2, lty = 2, color = 'blue'),
  
  pivot_longer(df, cols = everything()) %>% 
    ggplot(., aes(value, group = name)) +
    geom_histogram(aes(fill = name, color = name), position = position_identity(), alpha = .5)
)


ggplot(df, aes(ref_cbi)) +
  geom_density() +
  geom_boxplot(width = .2)

