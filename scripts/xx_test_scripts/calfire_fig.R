# figures for the calfire proposal

# concept figure for tree mortality on fire severity
# a) heat map of risk with trees overlaid
# b) distance vs. effect
# c) tri-variate surface plot with intensity + area vs. risk

rm(list = ls())

library(tidyverse)
library(sf)
library(terra)
library(cmdstanr)
library(bayesplot)
library(plotly)
library(cowplot)

# functions
HPDI <- rethinking::HPDI
zscore <- function(x) (x - mean(x)) / sd(x)
standarize <- function(x) x / max(x)
dens <- rethinking::dens

source('scripts/functions/functions_process_spatial.R')


# other functions ---------------------------------------------------------

f_plot_w <- function(wmatrix){
  plot(NA, xlim = range(dist_vec_sim), ylim = c(0, .05))
  for(i in 1:100){
    lines(dist_vec_sim, wmatrix[i,], col= scales::alpha(1, .2))
  }
}
f_wmax <- function(wmatrix){
  wmax <- dist_vec_sim[apply(w_rep, 1, which.max)]
  return(list(dens(wmax), HPDI(wmax, .9)))
}
f_w90 <- function(wmatrix){
  w90 <- map_vec(1:nrow(wmatrix), function(i){
    idx_90th <- which(cumsum(wmatrix[i,]) >= .90)[1]
    dist_vec_sim[idx_90th]
  })
  cumsum(wmatrix[1,]) %>% plot(dist_vec_sim, ., type = 'l', col= scales::alpha(1, .2),
                               xlab = 'Distance (km)', ylab = 'Cumulative effect')
  for(i in 1:50){
    cumsum(wmatrix[i,]) %>% lines(dist_vec_sim, ., col= scales::alpha(1, .2)) 
  }
  abline(h = .9, lty = 2, lwd = 2, col = 'red')
  segments(x0 = HPDI(w90, .9)[1], x1 = HPDI(w90, .9)[2], y0 = .9, y1 = .9,
           lwd = 4, col = 'red4')
  points(x = median(w90), y = .9, pch = 16, cex = 2, col = 'brown3')
  
  HPDI(w90, .9)
}

f_w90gg <- function(wmatrix){
  w90_df <- sapply(1:50, function(i) cumsum(wmatrix[i,])) %>% 
    as_tibble() %>% 
    mutate(distance = dist_vec_sim, .before = V1) %>% 
    pivot_longer(cols = starts_with('V'), names_prefix = 'V', names_to = 'draw')
  w90 <- map_vec(1:nrow(wmatrix), function(i){
    idx_90th <- which(cumsum(wmatrix[i,]) >= .90)[1]
    dist_vec_sim[idx_90th]
  })
  range_df <- tibble(value = .9, 
                     distance = median(w90), 
                     lower = HPDI(w90, .9)[1],
                     upper = HPDI(w90, .9)[2])
  
  ggplot(w90_df, aes(distance, value)) +
    geom_line(aes( group = draw), alpha = .1) +
    geom_hline(yintercept = .9, lty = 2, color = 'red3', lwd = .5, alpha = .5) +
    geom_pointrange(data = range_df, 
                    aes(distance, value, xmin = lower, xmax = upper), 
                    color = 'red4', size = .9, linewidth = 1)
  
}




# get ring data -----------------------------------------------------------

# get it from your actual data. pretend that ndvi is tree mortality.
dat <- read_rds('outputs/tabular/model_knp_ringdata.rds')


# make distance dependent data wide. N plots (rows) x M rings (cols)
# instaed of zscoring, make between 0 and 1. 
dd_var <- dat$rings %>%
  mutate(var_s = zscore(ndvi.mean)) %>% 
  pivot_wider(id_cols = id, names_from = radius, values_from = var_s, names_prefix = 'X_') %>% 
  arrange(id) %>% 
  select(-id) %>% 
  as.matrix()

# get area and distance vectors out. assume same for each focal site. 
dist_vec <- dat$rings %>% 
  mutate(radius = radius / 1000) %>% 
  filter(id == 1) %>% 
  pull(radius) 
inner <- c(0, dist_vec[-length(dist_vec)]) # inner ring radius
area_vec <- f_ring_area(dist_vec, inner) 



# simulate data -----------------------------------------------------------

set.seed(123)
# simulate response variable
#noise <- rlnorm(N)
delta <- .3
a0 <- 0
beta <- .5
sigma <- 1.1
N <- nrow(dd_var)
exp(-dist_vec^2 / (2 * delta^2)) %>% plot(dist_vec, ., ylim = c(0, 1))
w0 <- exp(-dist_vec^2 / (2 * delta^2)) * area_vec
w <- w0 / sum(w0) 
plot(dist_vec, w)
DD <- dd_var %*% w
mu <- a0 + beta * DD 
y <- rnorm(N, mu, sigma)

# generate distance and area vectors to simulate w to new data
lm_a_d <- lm(area_vec ~ dist_vec) # luckily there's a linear relationship
dist_vec_sim <- seq(.05, 1, by = .01)
area_vec_sim <- predict(lm_a_d, list(dist_vec = dist_vec_sim))





# create counterfactuals --------------------------------------------------


# choose one random focal site where y > 3
idx <- sample(which(y > 3), 1)



# apply counterfactual treatments on this site, adjusting dd_var
red_vec <- seq(0, 1, by = .05)
design <- expand_grid(reduction = red_vec,
                      distance = dist_vec) %>% 
  mutate(area = rep(cumsum(area_vec), times = length(red_vec)))


# get raw data before zscoring
dd_var_raw <- dat$rings %>%
  pivot_wider(id_cols = id, names_from = radius, values_from = ndvi.mean, names_prefix = 'X_') %>% 
  arrange(id) %>% 
  select(-id) %>% 
  as.matrix()

# get dd_var to simulate data from 
dd_var_0 <- dd_var_raw[idx,]
dd_var_l <- map2(design$distance, design$reduction, 
     ~ ifelse(dist_vec <= .x, dd_var_0 * (1-.y), dd_var_0))
dd_var_sim <- matrix(unlist(dd_var_l), byrow = T, 
                     nrow = nrow(design), ncol = length(dist_vec))

# zscore it now
dd_var_simz <- (dd_var_sim - mean(dat$rings$ndvi.mean)) / sd(dat$rings$ndvi.mean)






# model -------------------------------------------------------------------


# create data list for stan
dat_list <- list(
  N = length(y), 
  M = ncol(dd_var),
  dd_var = dd_var,
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y,
  Nsim_c = nrow(dd_var_simz),
  dd_var_sim = dd_var_simz,
  Nsim_d = length(dist_vec_sim),
  dist_sq_sim = dist_vec_sim^2,
  area_sim = area_vec_sim
)
str(dat_list)

# model with stan??!
curve(dlnorm(x, -1, 1.5), 0, 1) # playing with prior options for delta
curve(pscl::densigamma(x, 1, .2), .01, 1)


stanmod <- cmdstan_model('scripts/stan_models/basic_DD_Norm.stan')


fit01 <- stanmod$sample(data = dat_list, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 

fit01$summary(variables = c('delta', 'a0', 'beta', 'sigma')) %>% 
  print(n = Inf)
mcmc_pairs(fit01$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit01$draws(c( 'w[2]', 'w[3]', 'w[20]')))

ppc_dens_overlay(y, fit01$draws('y_rep', format = 'matrix')[1:50,])
w_fit <- fit01$draws('w', format = 'matrix')
ppc_intervals(w, yrep = w_fit)

# find wmax and w90
w_rep <- fit01$draws('w_rep', format = 'matrix')
f_plot_w(w_rep)
p_w90 <- f_w90gg(w_rep) +
  labs(x = 'Distance (km)', y = 'Cumulative effect of \nmortality on fire severity') +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(size = 10))
p_w90
f_wmax(w_rep)




# counterfactuals viz -----------------------------------------------------


# check out counterfactuals
y_sim <- fit01$draws('y_sim', format = 'matrix')
y_sim_risk <- rethinking::inv_logit(y_sim)

df_cf <- design %>% 
  mutate(mean = apply(y_sim, 2, mean),
         Risk = apply(y_sim_risk, 2, mean), 
         lower = apply(y_sim_risk, 2, HPDI, .9)[1,],
         upper = apply(y_sim_risk, 2, HPDI, .9)[2,], 
         ha = area * 100)
ggplot(df_cf, aes(distance, Risk, color = reduction)) + geom_point()

p <- plot_ly(df_cf, y = ~ha, x = ~reduction, z = ~ Risk, 
        type = 'mesh3d', intensity = ~Risk)  %>% hide_colorbar()
p
p %>% 
  layout(scene = list(xaxis = list(title = "Intensity", nticks = 2, titlefont = list(size = 45), nticks = 3, tickfont = list(size = 18)), 
                      yaxis = list(title = "Area treated", nticks = 1, titlefont = list(size = 45), nticks = 3, tickfont = list(size = 18)),        
                      zaxis = list(title = "Fire risk", nticks = 2, titlefont = list(size = 45), tickfont = list(size = 18))))
# saved to png interactively




# plot map ----------------------------------------------------------------

r_brick <- rast('outputs/spatial/compiled/rasters_knp.tif') 

# plot a focal site. choose a high severity spot
#set.seed(0)
#idx <- sample(which(y > 3), 1)

focalsite <- st_as_sf(dat$loc[idx,], coords = c('x', 'y'), crs = crs(r_brick))
point_buffer <- st_buffer(focalsite, 1250)
cbi_crop <- crop(r_brick$cbi, point_buffer)

# make rings around taht focal site
circles <- map(dist_vec*1000, ~ st_buffer(focalsite, .x)) 
rings <- bind_rows(
  circles[[1]],
  map_df(2:length(circles), function(i){
    suppressWarnings(
      st_difference(circles[[i]], circles[[i-1]]$geometry)
    )
  })
) %>% 
  mutate(radius = dist_vec*1000, .before = geometry, 
         alpha = radius / 1000) 


# make a mask
poly <- as.polygons(ext(cbi_crop)) %>% st_as_sf()
st_crs(poly) <- st_crs(rings)
outerring <- st_buffer(focalsite, 1000)
polymask <- st_difference(poly, outerring)


# make the map
rastdf <- spatSample(cbi_crop, 100000, xy = T)
p_map <- ggplot() +
  geom_tile(data = rastdf, aes(x, y, fill = cbi)) +
  geom_sf(data = polymask, color = NA, fill = 'grey20', alpha = .9) +
  geom_sf(data = rings, aes(alpha = alpha),fill = 'grey30', color = NA) +
  geom_sf(data = outerring, fill = NA, color = 'grey10', linewidth= 1) +
  geom_sf(data = focalsite, size = 4) +
  scale_alpha(range = c(0, .9)) +
  scale_fill_distiller(palette = 'RdYlGn') +
  theme_void() + 
  theme(legend.position = 'none', 
        panel.border = element_rect(colour = 'black', fill = NA, linewidth = 1)) 
p_map

# assemble everything -----------------------------------------------------
library(magick)



p_plotly <- ggdraw() + draw_image('figures/calfire_plotly.png')

p_title <- ggdraw() + 
  draw_label(fontface = 'bold', size = 18, fontfamily = "sans",
    "Hypothetical effect of tree mortality on fire severity"
  ) 


grid_R <- plot_grid(p_w90, p_plotly, nrow = 2,scale = c(.9, 1.1), 
          labels = c('B', 'C'), hjust = -1) 

grid_LR <- plot_grid(p_map, grid_R, ncol = 2, 
                        scale = c(.85, 1), labels = c('A',''), hjust = -1)

grid_final <- plot_grid(p_title, grid_LR, ncol = 1, rel_heights = c(.2, 1))


grid_final
ggsave('figures/calfire_final.png', grid_final, width = 9.4, height = 4.8)





