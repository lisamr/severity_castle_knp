library(tidyverse)
library(cmdstanr)
library(bayesplot)


dat <- read_rds('outputs/tabular/model_knp_data.rds')
HPDI <- rethinking::HPDI
zscore <- function(x) (x - mean(x)) / sd(x)
dens <- rethinking::dens



# make distance dependent data wide. N plots (rows) x M rings (cols)
dd_var <- dat$rings %>%
  filter(radius < 1000) %>% 
  mutate(var_z = zscore(sd)) %>% 
  pivot_wider(id_cols = id, names_from = radius, values_from = var_z, names_sep = '.') %>% 
  arrange(id) %>% 
  select(-id) %>% 
  as.matrix()
dd_var <- dd_var[-which(is.na(dat$preds$cbi)),]

# get area and distance vectors out. assume same for each focal site. 
area_vec <- dat$rings %>% 
  filter(id == 1) %>% 
  filter(radius < 1000) %>% 
  pull(ha) %>% 
  # might need to set back to m2 for consistent distance and area units
  units::set_units('km^2') %>% 
  as.numeric()
dist_vec <- dat$rings %>%
  filter(radius < 1000) %>% 
  mutate(radius = radius / 1000) %>% 
  filter(id == 1) %>% 
  pull(radius) 

# pull out response variable
cbi <- dat$preds$cbi
cbi <- cbi[!is.na(cbi)]
y <- as.integer(cbi > 1.5) # normally would be 2.25 but somethings wrong with the output
sum(y)

# generate distance and area vectors to simulate w to new data
lm_a_d <- lm(area_vec ~ dist_vec) # luckily there's a linear relationship
dist_vec_sim <- seq(.05, 2, by = .01)
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
stanmod <- cmdstan_model('scripts/stan_models/basic_DD_Bern.stan')

fit01 <- stanmod$sample(data = dat_list, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 

fit01$summary(variables = c('delta', 'a0', 'beta', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit01$draws(c('delta', 'a0', 'beta')))
mcmc_pairs(fit01$draws(c( 'w[2]', 'w[3]', 'w[4]')))
ppc_dens_overlay(y, fit01$draws('y_rep', format = 'matrix')[1:50,])



# find wmax and w90
w_rep <- fit01$draws('w_rep', format = 'matrix')
plot(NA, xlim = range(dist_vec_sim), ylim = c(0, .1))
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




# try with normal ---------------------------------------------------------

stanmod <- cmdstan_model('scripts/stan_models/basic_DD_Norm.stan')

dat_list_norm <- dat_list
dat_list_norm$y <- zscore(cbi)

fit02 <- stanmod$sample(data = dat_list_norm, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 
fit02$summary(variables = c('delta', 'a0', 'beta', 'sigma', 'w')) %>% 
  print(n = Inf)
mcmc_pairs(fit02$draws(c('delta', 'a0', 'beta')))
mcmc_pairs(fit02$draws(c( 'w[2]', 'w[3]', 'w[4]')))

# find wmax and w90
w_rep <- fit02$draws('w_rep', format = 'matrix')
plot(NA, xlim = range(dist_vec_sim), ylim = c(0, .1))
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

delta_draw <- fit02$draws('delta', format = 'matrix')
dens(delta_draw, lwd= 2)
x <- seq(.001, 12, by = .01)
dlnorm(x) %>% lines(x, ., lwd = 2, lty = 2)
