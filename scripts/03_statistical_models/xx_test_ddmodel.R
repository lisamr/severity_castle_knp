library(tidyverse)
library(cmdstanr)
library(bayesplot)

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
                               xlab = 'dist (km)', ylab = 'cumulative effect')
  for(i in 1:50){
    cumsum(wmatrix[i,]) %>% lines(dist_vec_sim, ., col= scales::alpha(1, .2)) 
  }
  abline(h = .9, lty = 2, lwd = 2, col = 'red')
  segments(x0 = HPDI(w90, .9)[1], x1 = HPDI(w90, .9)[2], y0 = .9, y1 = .9,
           lwd = 4, col = 'red4')
  points(x = median(w90), y = .9, pch = 16, cex = 2, col = 'brown3')
  
  HPDI(w90, .9)
}



# load --------------------------------------------------------------------


dat <- read_rds('outputs/tabular/model_knp_ringdata.rds')
HPDI <- rethinking::HPDI
zscore <- function(x) (x - mean(x)) / sd(x)
dens <- rethinking::dens
dat$rings %>% head


# make distance dependent data wide. N plots (rows) x M rings (cols)
dd_var <- dat$rings %>%
  mutate(var_z = zscore(ndvi.sd)) %>% 
  pivot_wider(id_cols = id, names_from = radius, values_from = var_z, names_prefix = 'X_') %>% 
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


# simulate response variable
#noise <- rlnorm(N)
delta <- .2
a0 <- -.2
beta <- 1
sigma <- .5
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

# create data list for stan
dat_list <- list(
  N = length(y), 
  M = ncol(dd_var),
  dd_var = dd_var,
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y,
  Nsim = length(dist_vec_sim),
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
f_w90(w_rep)
f_wmax(w_rep)

# try with real data ------------------------------------------------------


rmvals <- is.na(dat$preds$cbi)

dd_var2 <- dd_var[!rmvals,]
y2 <- zscore(dat$preds$cbi[!rmvals])

dat_list2 <- list(
  N = length(y2), 
  M = ncol(dd_var2),
  dd_var = dd_var2,
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y2,
  Nsim = length(dist_vec_sim),
  dist_sq_sim = dist_vec_sim^2,
  area_sim = area_vec_sim
)
str(dat_list2)

fit02 <- stanmod$sample(data = dat_list2, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 
fit02$summary(variables = c('delta', 'a0', 'beta', 'sigma')) %>% 
  print(n = Inf)
mcmc_pairs(fit02$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit02$draws(c( 'w[2]', 'w[3]', 'w[20]')))


# find wmax and w90
w_rep <- fit02$draws('w_rep', format = 'matrix')
f_plot_w(w_rep)
f_w90(w_rep)
f_wmax(w_rep)

# furhter sample to see how sample size effects things
keepvals <- sample(which(!rmvals), 100)
dd_var3 <- dd_var[keepvals,]
y3 <- zscore(dat$preds$mtbs_CBI[keepvals])

dat_list3 <- list(
  N = length(y3), 
  M = ncol(dd_var3),
  dd_var = dd_var3,
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y3,
  Nsim = length(dist_vec_sim),
  dist_sq_sim = dist_vec_sim^2,
  area_sim = area_vec_sim
)


fit03 <- stanmod$sample(data = dat_list3, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 
fit03$summary(variables = c('delta', 'a0', 'beta', 'sigma', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit03$draws(c('delta', 'a0', 'beta', 'sigma')))

# find wmax and w90
w_rep <- fit03$draws('w_rep', format = 'matrix')
f_plot_w(w_rep)
f_w90(w_rep)
f_wmax(w_rep)


cowplot::plot_grid(
  mcmc_pairs(fit03$draws(c( 'w[2]', 'w[3]', 'w[20]'))),
  mcmc_pairs(fit02$draws(c( 'w[2]', 'w[3]', 'w[20]'))),
  labels = c('n = 100 focal sites', 'n = 1000'), scale = .9
)


