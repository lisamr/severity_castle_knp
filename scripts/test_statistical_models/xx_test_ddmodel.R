library(tidyverse)
library(cmdstanr)
library(bayesplot)

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_distance_dependence.R')


# load --------------------------------------------------------------------


dat <- read_rds('outputs/tabular/model_knp_ringdata.rds')
HPDI <- rethinking::HPDI
zscore <- function(x) (x - mean(x)) / sd(x)
dens <- rethinking::dens
dat$rings %>% head
# make sure no NAs
apply(dat$rings, 2, function(x) sum(is.na(x)))


# make a list of distance dependent data wide. N plots (rows) x M rings (cols) X Nlayers
variables <- names(dat$rings)[-(1:2)]
ddvar_wide_list <- lapply(variables, f_ddvar_wide) %>% set_names(., variables)




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
N <- nrow(ddvar_wide_list[[1]])
exp(-dist_vec^2 / (2 * delta^2)) %>% plot(dist_vec, ., ylim = c(0, 1))
w0 <- exp(-dist_vec^2 / (2 * delta^2)) * area_vec
w <- w0 / sum(w0) 
plot(dist_vec, w)
DD <- ddvar_wide_list[[1]] %*% w
mu <- a0 + beta * DD 
y <- rnorm(N, mu, sigma)

# generate distance and area vectors to simulate w to new data
lm_a_d <- lm(area_vec ~ dist_vec) # luckily there's a linear relationship
dist_vec_sim <- seq(.05, 1, by = .01)
area_vec_sim <- predict(lm_a_d, list(dist_vec = dist_vec_sim))

# create data list for stan
dat_list <- list(
  N = length(y), 
  M = ncol(ddvar_wide_list[[1]]),
  dd_var = ddvar_wide_list[[1]],
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y,
  # for simulating counterfactuals and distance functions
  Nsim_c = 1L, 
  dd_var_sim = rbind(rep(0, ncol(ddvar_wide_list[[1]]))), 
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
#mcmc_pairs(fit01$draws(c( 'w[2]', 'w[3]', 'w[20]')))

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

dd_var2 <- ddvar_wide_list[[3]][!rmvals,]
y2 <- zscore(dat$preds$cbi[!rmvals])

dat_list2 <- list(
  N = length(y2), 
  M = ncol(dd_var2),
  dd_var = dd_var2,
  area = area_vec,
  dist_sq = dist_vec^2,
  y = y2,
  # for simulating counterfactuals and distance functions
  Nsim_c = 1L, 
  dd_var_sim = rbind(rep(0, ncol(dd_var2))), 
  Nsim_d = length(dist_vec_sim),
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
ppc_dens_overlay(y2, fit02$draws('y_rep', format = 'matrix')[1:50,])

# find wmax and w90
w_rep <- fit02$draws('w_rep', format = 'matrix')
f_plot_w(w_rep)
f_w90(w_rep)
f_wmax(w_rep)



# why are you getting errors ----------------------------------------------

# sampler is complaining when w0=0. leads to w being NaN. 
f_decay <- function(d_sq, a, delt) exp(-d_sq / (2 * (delt^2))) * a
tmp <- f_decay(dist_vec_sim^2, area_vec_sim, .03)
tmp[1]
sum(tmp)
plot(dist_vec_sim, tmp, type = 'l')

# try out bernoulli likelihood --------------------------------------------

stanmod_bern <- cmdstan_model('scripts/stan_models/basic_DD_Bern.stan')

dat_list3 <- dat_list2
dat_list3$y <- as.integer(dat$preds$cbi[!rmvals] >= 2.25)

fit03 <- stanmod_bern$sample(data = dat_list3, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4,
                        chains = 4) 

fit03$summary(variables = c('delta', 'a0', 'beta')) %>% 
  print(n = Inf)
fit03$summary(variables = c('w', 'w0')) %>% print(n = Inf)

mcmc_pairs(fit03$draws(c('delta', 'a0', 'beta')))
mcmc_pairs(fit03$draws(c( 'w[2]', 'w[3]', 'w[20]')))

ppc_dens_overlay(dat_list3$y, fit03$draws('y_rep', format = 'matrix')[1:50,])

# find wmax and w90
w_rep <- fit03$draws('w_rep', format = 'matrix')
f_plot_w(w_rep)
f_w90(w_rep)
f_wmax(w_rep)





# limit sample size -------------------------------------------------------






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


