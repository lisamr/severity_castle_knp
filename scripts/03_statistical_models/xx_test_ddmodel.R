library(tidyverse)
library(cmdstanr)
library(bayesplot)

source('scripts/functions/functions_process_spatial.R')

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
stanmod <- cmdstan_model('scripts/stan_models/basic_DD_Norm.stan')

fit01 <- stanmod$sample(data = dat_list, 
               iter_warmup = 1000,
               iter_sampling = 1000,
               parallel_chains = 4) 

fit01$summary(variables = c('delta', 'a0', 'beta', 'sigma', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit01$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit01$draws(c( 'w[2]', 'w[3]', 'w[20]')))

ppc_dens_overlay(y, fit01$draws('y_rep', format = 'matrix')[1:50,])
w_fit <- fit01$draws('w', format = 'matrix')
ppc_intervals(w, yrep = w_fit)

# find wmax and w90
w_rep <- fit01$draws('w_rep', format = 'matrix')
plot(NA, xlim = range(dist_vec_sim), ylim = c(0, .05))
for(i in 1:100){
  lines(dist_vec_sim, w_rep[i,], col= scales::alpha(1, .2))
}
wmax <- dist_vec_sim[apply(w_rep, 1, which.max)]
dens(wmax); HPDI(wmax, .9)
w90 <- map_vec(1:nrow(w_rep), function(i){
  idx_90th <- which(cumsum(w_rep[i,]) >= .90)[1]
  dist_vec_sim[idx_90th]
})
dens(w90); HPDI(w90, .9)
cumsum(w_rep[1,]) %>% plot(type = 'l')


# try with real data ------------------------------------------------------


rmvals <- is.na(dat$preds$mtbs_CBI)

dd_var2 <- dd_var[!rmvals,]
y2 <- zscore(dat$preds$mtbs_CBI[!rmvals])

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


fit02 <- stanmod$sample(data = dat_list2, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 
fit02$summary(variables = c('delta', 'a0', 'beta', 'sigma', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit01$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit02$draws(c( 'w[2]', 'w[3]', 'w[20]')))


# find wmax and w90
w_rep <- fit02$draws('w_rep', format = 'matrix')
plot(NA, xlim = range(dist_vec_sim), ylim = c(0, .05))
for(i in 1:100){
  lines(dist_vec_sim, w_rep[i,], col= scales::alpha(1, .2))
}
wmax <- dist_vec_sim[apply(w_rep, 1, which.max)]
dens(wmax); HPDI(wmax, .9)
w90 <- map_vec(1:nrow(w_rep), function(i){
  idx_90th <- which(cumsum(w_rep[i,]) >= .90)[1]
  dist_vec_sim[idx_90th]
})
dens(w90); HPDI(w90, .9)
cumsum(w_rep[1,]) %>% plot(dist_vec_sim, ., type = 'l', col= scales::alpha(1, .2),
                           xlab = 'dist (km)', ylab = 'cumulative effect')
for(i in 1:50){
  cumsum(w_rep[i,]) %>% lines(dist_vec_sim, ., col= scales::alpha(1, .2)) 
}
abline(h = .9, lty = 2, lwd = 2, col = 'red')
