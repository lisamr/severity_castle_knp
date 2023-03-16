# cbi data from GEE seems really off. ~3% of knp > 2.25. 
# compare to mtbs

library(terra)
library(tidyverse)


# load data
cbi <- rast('outputs/spatial/GEE_CBI/issues/KNP_2021_CBI_bc.tif')
cbi_CT <- rast('outputs/spatial/GEE_CBI/issues/KNP_2021_CBI_bc_CT.tif')
cbi_l9 <- rast('outputs/spatial/GEE_CBI/issues/KNP_2021_CBI_bc_l9.tif')
cbi_cp <- rast('outputs/spatial/GEE_CBI/issues/KNP_2021_CBI_bc_copy.tif')
cbi_cp2 <- rast('outputs/spatial/GEE_CBI/issues/KNP_2021_CBI_bc_copy2.tif')
cbi_cor <- rast('outputs/spatial/GEE_CBI/KNP_2021_CBI_bc_corrected.tif')
mtbs <- rast('data/fires/mtbs_SQF_KNP.grd')



#this function was used to convert rdnbr to cbi and ba
# extendedAssessments <- function(RdNBR){
#   RdNBR[RdNBR <= -369] <- -368.99999999
#   CBI <- 1/.388 * log((RdNBR + 369)/421.7)
#   CBI[CBI < 0] <- 0
#   CBI[CBI > 3] <- 3
# 
#   # need to deal with extreme values.
#   RdNBR[RdNBR < 166.5] <- 166.5
#   RdNBR[RdNBR > 389*2] <- 389*2
#   BA <- (sin((RdNBR - 166.5)/389))^2 * 100
#   br <- rast(list(RdNBR, CBI, BA))
#   names(br) <- c('RdNBR', 'CBI', 'BA')
#   return(br)
# }

# force rasters to common extent. use gee layer.
mtbs_c <- crop(mtbs, cbi) 
mtbs_c$CBI[mtbs_c$CBI > 3] <- 3
brick <- list(cbi, CBI_CT =  cbi_CT, CBI_l9 = cbi_l9, CBI_cp = cbi_cp, 
              CBI_cp2 = cbi_cp2, CBI_cor = cbi_cor,
              mtbs = mtbs_c$CBI) %>% rast

# compare
brick_pts <- spatSample(brick, 5000, na.rm = T)
plot(brick_pts$CBI_bc, brick_pts$mtbs,  xlab = 'GEE cbi', ylab = 'ref cbi', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot( brick_pts$mtbs, brick_pts$CBI_CT,xlab = 'ref cbi',   ylab = 'GEE cbi CT', 
      pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot( brick_pts$mtbs, brick_pts$CBI_cp2, xlab = 'ref cbi',   ylab = 'GEE cbi copy', 
      pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot( brick_pts$mtbs, brick_pts$CBI_cor,xlab = 'ref cbi',   ylab = 'GEE cbi scale cor', 
      pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot(brick_pts$CBI_bc, brick_pts$CBI_CT, xlab = 'GEE cbi', ylab = 'GEE cbi CT', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot(brick_pts$CBI_bc, brick_pts$CBI_l9, xlab = 'GEE cbi', ylab = 'GEE cbi L9', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot(brick_pts$CBI_l9, brick_pts$CBI_cp, xlab = 'GEE cbi L9', ylab = 'GEE cbi copy', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot(brick_pts$CBI_bc, brick_pts$CBI_cp, xlab = 'GEE cbi', ylab = 'GEE cbi copy', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot(brick_pts$CBI_l9, brick_pts$CBI_cp2, xlab = 'GEE cbi L9', ylab = 'GEE cbi copy2', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)

plot(brick_pts$CBI_CT, brick_pts$CBI_cp, xlab = 'GEE cbi CT', ylab = 'GEE cbi copy', 
     pch = 16, col = scales::alpha(1, .2))
abline(0, 1, col = 'red', lty = 2, lwd = 2)




ggplot(brick_pts, aes(mtbs, CBI_cor)) +
  geom_point(alpha = .5) +
  geom_abline(slope = 1, intercept = 0, lwd = 2, color = 'lightblue') +
  coord_equal()


name_labels <- c(
  CBI_bc = 'L4-8, C02/T1_L2',
  CBI_l9 = 'L4-9, C02/T1_L2',
  CBI_cor = 'L4-9, C02/T1_L2\n corrected scale',
  CBI_cp = 'L4-7 C01/T1_SR,\n L8-9, C02/T1_L2',
  CBI_CT = 'L4-8 C01/T1_SR,\n L8-9, C02/T1_L2',
  mtbs = 'mtbs'
)
pivot_longer(brick_pts, cols = everything()) %>% 
  filter(!name %in% c('comp', 'CBI_cp2')) %>% 
  mutate(name = fct_relevel(name,names(name_labels) )) %>% 
  ggplot(., aes(value)) +
  geom_histogram(binwidth = .1) +
  facet_grid(rows = vars(name), labeller = labeller(name = name_labels))

ggplot(brick_pts, aes(CBI_CT)) +
  geom_density() +
  geom_boxplot(width = .2)


2.75e-05
0.0001	
