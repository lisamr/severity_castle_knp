library(tidyverse)
library(cmdstanr)
library(bayesplot)


dat <- read_rds('outputs/tabular/model_knp_data.rds')
HPDI <- rethinking::HPDI
zscore <- function(x) (x - mean(x)) / sd(x)
dens <- rethinking::dens
dat$rings


# make distance dependent data wide. N plots (rows) x M rings (cols)
dd_var <- dat$rings %>%
  mutate(means_z = zscore(mean)) %>% 
  pivot_wider(id_cols = id, names_from = radius, values_from = means_z, names_prefix = 'X_') %>% 
  arrange(id) %>% 
  select(-id) %>% 
  as.matrix()

# get area and distance vectors out. assume same for each focal site. 
area_vec <- dat$rings %>% 
  filter(id == 1) %>% 
  pull(ha) %>% 
  # might need to set back to m2 for consistent distance and area units
  units::set_units('km^2') %>% 
  as.numeric()
dist_vec <- dat$rings %>% 
  mutate(radius = radius / 1000) %>% 
  filter(id == 1) %>% 
  pull(radius) 


# simulate response variable
noise <- rlnorm(N)
delta <- 2 + noise # length scale ~ 500 m
a0 <- -.2
beta <- 1
sigma <- .5
N <- nrow(dd_var)
exp(-dist_vec^2 / (2 * 1^2)) %>% plot(ylim = c(0, 1))
w0 <- exp(-dist_vec^2 / (2 * delta^2)) * area_vec
w <- w0 / sum(w0) 
DD <- dd_var %*% w
mu <- a0 + beta * DD 
y <- rnorm(N, mu, sigma)

# generate distance and area vectors to simulate w to new data
lm_a_d <- lm(area_vec ~ dist_vec) # luckily there's a linear relationship
dist_vec_sim <- seq(.05, 2, by = .1)
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
plot(NA, xlim = range(dist_vec_sim), ylim = c(0, .02))
for(i in 1:50){
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

# fix first intercept at 0 ------------------------------------------------



# simulate response variable
delta <- .5 # length scale ~ 500 m
a0 <- -.2
beta <- 1
sigma <- .5
N <- nrow(dd_var)
dist_vec0 <- dist_vec[-1]
area_vec0 <- area_vec[-1]
w0 <- exp(-dist_vec0^2 / (2 * delta^2)) * area_vec0
w0 <- w0 / sum(w0) 
w <- c(0, w0)
DD <- dd_var %*% w
mu <- a0 + beta * DD
y <- rnorm(N, mu, sigma)

# create data list for stan
dat_list <- list(
  N = length(y), 
  M = ncol(dd_var),
  dd_var = dd_var,
  area = area_vec0,
  dist_sq = dist_vec0^2,
  y = y
)
str(dat_list)

# model with stan??!
stanmod <- cmdstan_model('scripts/stan_models/basic_DD_Norm_fix0.stan')

fit01.1 <- stanmod$sample(data = dat_list, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 

fit01.1$summary(variables = c('delta', 'a0', 'beta', 'sigma', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit01.1$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit01.1$draws(c( 'w[2]', 'w[3]', 'w[20]')))
mcmc_pairs(fit01.1$draws(c( 'w0[2]', 'w0[3]', 'w0[19]')))




# another try -------------------------------------------------------------

# what about nonspatial model? can you get that to work?

# create data list for stan
dat_list_2 <- list(
  N = length(y), 
  dd_var = dd_var[,1],
  y = y
)

# model with stan??!
stanmod_2 <- cmdstan_model('scripts/stan_models/nonspatial.stan')

fit02 <- stanmod_2$sample(data = dat_list_2, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4, seed = 1) 

fit02$summary(c('a0', 'beta','sigma'))
mcmc_pairs(fit02$draws(c('a0', 'beta','sigma')))


# try with simple simplex? ------------------------------------------------

# model with stan??!
stanmod_3 <- cmdstan_model('scripts/stan_models/simplex.stan')

fit03 <- stanmod_3$sample(data = dat_list, 
                          iter_warmup = 1000,
                          iter_sampling = 1000,
                          parallel_chains = 4) 

fit03$summary(variables = c('a0', 'beta', 'sigma', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit03$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit03$draws(c('w[1]', 'w[2]', 'w[3]', 'w[4]')))

ppc_dens_overlay(y, fit01$draws('y_rep', format = 'matrix')[1:50,])


