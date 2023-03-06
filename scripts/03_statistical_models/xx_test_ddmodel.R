library(tidyverse)
library(cmdstanr)
library(bayesplot)


dat <- read_rds('outputs/tabular/model_knp_data.rds')

zscore <- function(x) (x - mean(x)) / sd(x)
dens <- rethinking::dens
dat$rings
dat$preds %>% head

y <- dat$preds$cbi %>% zscore

# make distance dependent data wide. N plots (rows) x M rings (cols)
dd_var <- dat$rings %>%
  mutate(means_z = zscore(means)) %>% 
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



# create data list for stan
dat_list <- list(
  N = length(y), 
  M = ncol(dd_var),
  dd_var = dd_var,
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y
)
str(dat_list)

# model with stan??!
stanmod <- cmdstan_model('scripts/stan_models/basic_DD_Norm.stan')

fit01 <- stanmod$sample(data = dat_list, 
               iter_warmup = 1000,
               iter_sampling = 1000,
               parallel_chains = 4) 

fit01$summary(variables = c('delta', 'beta', 'sigma', 'w0')) %>% 
  print(n = Inf)
mcmc_pairs(fit01$draws(c('delta', 'a0', 'beta', 'sigma')))
mcmc_pairs(fit01$draws(c( 'w[2]', 'w[3]', 'w[20]')))
mcmc_pairs(fit01$draws(c( 'w0[2]', 'w0[3]', 'w0[20]')))

ppc_dens_overlay(y, fit01$draws('y_rep', format = 'matrix')[1:50,])



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


