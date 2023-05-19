# test out ordinal regressions on the CBI data

library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(tidybayes)

rm(list = ls())

dat <- read_rds('outputs/tabular/model_knp_ringdata.rds')
HPDI <- rethinking::HPDI
inv_logit <- rethinking::inv_logit
zscore <- function(x) (x - mean(x)) / sd(x)
zscore_sim <- function(xsim, x) (xsim - mean(x)) / sd(x) # for zscoring simulated data

# look at data
dat$preds %>% head
dat$preds$cbi %>% hist
summary(dat$preds$cbi)

# assign an ordered, categorical variable for CBI
# cutpoints from Miller & Thode 2007 (unchanged, low, mod, high)
cuts <- c(0, .1, 1.25, 2.25, 3)
dat$preds$cbi_factor <- cut(dat$preds$cbi, cuts, 
                labels = F, right = F, include.lowest = T)
dat$preds %>% 
  ggplot(aes(x = cbi_factor, y = cbi, group = cbi_factor)) +
  geom_boxplot()
dat$preds %>% 
  ggplot(aes(cbi, fill = as.factor(cbi_factor))) +
  geom_histogram(binwidth = .01) +
  harrypotter::scale_fill_hp_d(option = 'Ravenclaw')


prop_cbi <- table(dat$preds$cbi_factor) / sum(!is.na(dat$preds$cbi_factor))
cumsum(prop_cbi)





# add in simple predictors ------------------------------------------------

# load up the model
stanmod <- cmdstan_model('scripts/stan_models/ordinal.stan')

# use dead_area and ndvi as a predictors? rescale values too.
head(dat$preds)
dat1 <- dat$preds %>% 
  mutate(dead_area = dead_area/900, # %dead area now
         ndvi = ndvi/max(ndvi)
         )

# select data
dat2 <- dat1 %>% 
  as_tibble() %>% 
  select(cbi_factor, dead_area, ndvi) %>% 
  drop_na() %>% 
  mutate(across(-cbi_factor, zscore)) %>% 
  sample_frac(.25)
dat2X <- as.matrix(dat2[,-1])

# create data to simulate from 
dead_area_sim <- seq(0, 1, by = .1)
dat_sim <- tibble(dead_area = zscore_sim(dead_area_sim, dat1$dead_area),
       ndvi = 0) %>% as.matrix

dat_list <- list(
  N = nrow(dat2),
  K = max(dat2$cbi_factor),
  J = ncol(dat2X),
  X = dat2X,
  y = dat2$cbi_factor, 
  Nsim = nrow(dat_sim),
  Xsim = dat_sim
)
str(dat_list)


# fit the model
fit1 <- stanmod$sample(dat_list, parallel_chains = 4, chains = 4)



# check out results
fit1$summary(variables = c('B','c'))
# variable   mean median     sd    mad     q5    q95  rhat ess_bulk ess_tail
# <chr>     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
# 1 B[1]      0.232  0.232 0.0394 0.0388  0.168  0.297  1.00    3672.    2663.
# 2 B[2]      0.296  0.296 0.0399 0.0403  0.232  0.362  1.00    4198.    2782.
# 3 c[1]     -3.24  -3.24  0.102  0.100  -3.41  -3.07   1.00    2524.    2198.
# 4 c[2]     -1.11  -1.11  0.0468 0.0471 -1.19  -1.03   1.00    4150.    3238.
# 5 c[3]      0.666  0.667 0.0423 0.0416  0.596  0.733  1.00    4442.    3065.
draws <- fit1$draws(variables = c('B', 'c'), format = 'matrix')
draws_c <- fit1$draws('c', format = 'matrix')
inv_logit(apply(draws_c, 2, median))

# look at fit
y_rep <- fit1$draws('y_rep', format = 'matrix')
ppc_dens_overlay(dat_list$y, y_rep[1:50,])
ppc_bars(dat_list$y, y_rep[1:50,])
mcmc_pairs(fit1$draws(c('B', 'c')))

# make sure the intercepts make sense
cbi_counts <- table(dat2$cbi_factor)
prop_cbi <- cbi_counts/sum(cbi_counts)
as_tibble(inv_logit(draws_c)) %>%
  rename_with(~ as.character(parse_number(.x)), everything()) %>% 
  pivot_longer(everything()) %>% 
  ggplot(aes(value, group = name)) +
  stat_halfeye() +
  geom_vline(xintercept = cumsum(prop_cbi)[1:3], 
             color = 'slateblue', lty = 2, lwd = 1)
  
# get proportion of each CBI class? 
cum_prop <- inv_logit(draws_c) %>% cbind(0, ., 1)
post_prop <- apply(cum_prop, 1, diff) %>% 
  t %>% as_tibble %>% 
  setNames(1:4) %>% 
  pivot_longer(everything())
ggplot(post_prop, aes(value, fill = name, color = name)) +
  geom_histogram(bins = 100, alpha= .5) +
  geom_vline(xintercept = (prop_cbi), 
             color = 'slateblue', lty = 2, lwd = 1)



# vizualize counterfactuals -----------------------------------------------


# how to make sense of the predictors?? 
# I want to calculate the proportion of each class under different scenarios. 
phi_sim <- fit1$draws('phi_sim', format = 'matrix')

# cumulative proportions, in list form. 
x <- lapply(1:ncol(draws_c), function(k){
  # equivalent to: inv_logit(c - phi)
  inv_logit(sweep(-phi_sim, 1, draws_c[,k], FUN = '+')) 
})


# get posterior proportions for simulated data. 
post_sim <- lapply(seq_len(ncol(x[[1]])), 
                   # combine the columns of each matrix across the list
                   function(i) do.call(cbind, lapply(x, '[', , i))) %>% 
  # convert to proportions in each CBI class
  map(~ cbind(0, .x, 1)) %>% 
  map(~ apply(.x, 1, diff) %>% t) %>% 
  map(~ tibble(class = 1:(ncol(.x)),
               median = apply(.x, 2, median),
               lower = apply(.x, 2, HPDI, .9)[1,],
               upper = apply(.x, 2, HPDI, .9)[2,])) %>% 
  map2(.x = ., .y = 1:11, ~ cbind(.x, index = .y)) %>% 
  bind_rows() %>% 
  # merge with dat_sim
  left_join(as_tibble(dat_sim) %>% 
            mutate(dead_area = dead_area_sim,index = 1:nrow(.))
          )
ggplot(post_sim, aes(dead_area, median, color = class, fill = class, group = class)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA,alpha = .4) +
  labs(x = '% dead canopy', y = 'Probability of each class')




# add in a distance dependent variable ------------------------------------

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_distance_dependence.R')



# prep dist-dep variable. zscored in `f_ddvar_wide`
NA_cbi <- which(is.na(dat$preds$cbi))
dat$rings <- dat$rings %>% 
  mutate(dead_area.mean = dead_area.mean/900) %>% # omit once you fix this
  filter(!(id %in% NA_cbi))
ddvar_wide <- f_ddvar_wide(dat$rings, 'dead_area.mean') # dead_area.mean

# get area and distance vectors out
dist_vec <- dat$rings %>% 
  mutate(radius = radius / 1000) %>% 
  filter(id == 1) %>% 
  pull(radius) 
inner <- c(0, dist_vec[-length(dist_vec)]) # inner ring radius
area_vec <- f_ring_area(dist_vec, inner) 

# other fixed effects
dat3 <- dat$preds %>% 
  filter(!(ID %in% NA_cbi)) %>% 
  select(cbi_factor, ndvi, elevation) %>% 
  mutate(across(-cbi_factor, zscore))
fixed_effects <- select(dat3, -cbi_factor) %>% as.matrix

# create data to simulate from 
dat_sim <- tibble(ndvi = 0) %>% as.matrix

dat_list <- list(
  N = nrow(dat3),
  K = max(dat3$cbi_factor),
  J = ncol(fixed_effects),
  X = fixed_effects,
  y = dat3$cbi_factor,
  M = ncol(ddvar_wide),
  dd_var = ddvar_wide,
  area = area_vec,
  dist_sq = dist_vec
  # Nsim = nrow(dat_sim),
  # Xsim = dat_sim,
  # dd_var_sim = 
)
str(dat_list)


# load up the model
stanmod_DD <- cmdstan_model('scripts/stan_models/ordinal_DD.stan')

# fit the model
fit_dd <- stanmod_DD$sample(dat_list, parallel_chains = 4, chains = 4, 
                            iter_sampling = 1000, iter_warmup = 1000)
fit_dd$summary(variables = c('B', 'beta', 'c', 'delta'))
# # A tibble: 6 × 10
# variable   mean median     sd    mad     q5    q95  rhat ess_bulk ess_tail
# <chr>     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
# 1 B[1]      0.295  0.294 0.0194 0.0192  0.262  0.326 1.00     4987.    2701.
# 2 beta      0.342  0.342 0.0187 0.0189  0.311  0.373 1.00     3781.    2994.
# 3 c[1]     -3.27  -3.27  0.0505 0.0502 -3.35  -3.18  1.00     2842.    2467.
# 4 c[2]     -1.07  -1.07  0.0227 0.0227 -1.11  -1.03  0.999    4366.    3620.
# 5 c[3]      0.714  0.714 0.0218 0.0214  0.678  0.749 1.00     3892.    3068.
# 6 delta     0.129  0.129 0.0116 0.0111  0.111  0.149 1.00     4379.    2804.

# look at fit
y_rep <- fit_dd$draws('y_rep', format = 'matrix')
ppc_dens_overlay(dat_list$y, y_rep[1:50,])
ppc_bars(dat_list$y, y_rep[1:50,])
mcmc_pairs(fit_dd$draws(c('B', 'c', 'beta', 'delta')))



# include multiple DD variables -------------------------------------------

# i adapted the stan model to include an array of matrices
ddvar_wide_array <- map(c('count_red.mean', 'count_grey.mean'), 
                        ~f_ddvar_wide(dat$rings, .x))


dat_list <- list(
  N = nrow(dat3),
  K = max(dat3$cbi_factor),
  J = ncol(fixed_effects),
  X = fixed_effects,
  y = dat3$cbi_factor,
  nDD = length(ddvar_wide_array),
  M = ncol(ddvar_wide),
  dd_var = ddvar_wide_array,
  area = area_vec,
  dist_sq = dist_vec
  # Nsim = nrow(dat_sim),
  # Xsim = dat_sim,
  # dd_var_sim = 
)
str(dat_list)

# load up the model
stanmod_DDarray <- cmdstan_model('scripts/stan_models/ordinal_DDarray.stan')

# fit the model
fit_ddarray <- stanmod_DDarray$sample(dat_list, parallel_chains = 4, chains = 4, 
                            iter_sampling = 1000, iter_warmup = 1000)
fit_ddarray$summary(variables = c('B', 'beta', 'c', 'delta'))



# 
# # model it ----------------------------------------------------------------
# 
# # use an ordinal regression to estimate the intercepts of the cutpoints. 
# # will give you proportions of each class. 
# 
# stanmod <- cmdstan_model('scripts/stan_models/ordinal.stan')
# 
# # remove NAs, sample 25% of the data so it goes faster
# cbi_factor2 <- dat$preds$cbi_factor[!is.na(dat$preds$cbi_factor)]
# cbi_factor2 <- sample(cbi_factor2, length(cbi_factor2)/4)
# 
# dat_list <- list(
#   N = length(cbi_factor2),
#   K = max(cbi_factor2),
#   y = cbi_factor2
# )
# str(dat_list)
# 
# # model it
# fit1 <- stanmod$sample(dat_list, parallel_chains = 4, chains = 4)
# 
# # check out results
# fit1$summary(variables = c('c'))
# draws <- fit1$draws(variables = c('c'), format = 'matrix')
# inv_logit(apply(draws, 2, median))
# 
# # look at fit
# y_rep <- fit1$draws('y_ppc', format = 'matrix')
# ppc_dens_overlay(cbi_factor2, y_rep[1:50,])
# ppc_bars(cbi_factor2, y_rep[1:50,])
# 
# datprop <- tibble(name = 1:3, value = cumsum(prop_cbi)[1:3])
# as_tibble(draws) %>% 
#   mutate(gamma = 0,
#          across(`c[1]`:`c[3]`, ~ inv_logit(.x - gamma)),
#          gamma = NULL) %>% 
#   rename_with(~ as.character(parse_number(.x)), everything()) %>% 
#   pivot_longer(everything()) %>% 
#   ggplot(aes(name, value)) +
#   stat_pointinterval() +
#   geom_point(data = datprop, color = 'slateblue1', pch = '*', size = 10)
# 

