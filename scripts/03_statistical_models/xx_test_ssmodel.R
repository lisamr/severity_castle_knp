rm(list = ls())

library(tidyverse)
library(cmdstanr)
library(bayesplot)


dat <- read_rds('outputs/tabular/model_knp_bufferdata.rds')
HPDI <- rethinking::HPDI
zscore <- function(x) (x - mean(x)) / sd(x)
dens <- rethinking::dens


# make distance dependent data wide. N plots (rows) x M rings (cols)
dd_var <- dat$buffers %>%
  mutate(means_z = zscore(mean)) %>% 
  pivot_wider(id_cols = id, names_from = radius, values_from = means_z, names_prefix = 'X_') %>% 
  arrange(id) %>% 
  select(-id) %>% 
  as.matrix()


# simulate response variable
ss <- 3.4
a0 <- -.2
beta <- 1
sigma <- .5
N <- nrow(dd_var)
# linearly interpolate var
f_interpol <- function(SS, E) {
  if (SS > ncol(E)) {
    return(E[, ncol(E)])
  }
  SS_trunc = as.integer(trunc(SS))
  prop_scales = SS - SS_trunc
  E[,SS_trunc] * (1 - prop_scales) + E[,SS_trunc + 1] * prop_scales
}
w <- f_interpol(ss, dd_var)
# put it together
mu <- a0 + beta * w 
y <- rnorm(N, mu, sigma)


# create data list for stan
dat_list <- list(
  N = length(y), 
  M = ncol(dd_var),
  dd_var = dd_var,
  y = y
)
str(dat_list)

# model with stan??!
stanmod <- cmdstan_model('scripts/stan_models/basic_ss_Norm.stan')

fit01 <- stanmod$sample(data = dat_list, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 

fit01$summary(variables = c('ss', 'a0', 'beta', 'sigma')) %>% 
  print(n = Inf)
mcmc_pairs(fit01$draws(c('ss', 'a0', 'beta', 'sigma')))


# try with real data? -----------------------------------------------------

y <- zscore(dat$preds$cbi)
dat_list2 <- dat_list
dat_list2$y <- y

fit02 <- stanmod$sample(data = dat_list2, 
                        iter_warmup = 1000,
                        iter_sampling = 1000,
                        parallel_chains = 4) 
fit02$summary(variables = c('ss', 'a0', 'beta', 'sigma'))
mcmc_pairs(fit02$draws(c('ss', 'a0', 'beta', 'sigma')))

curve(dlnorm(x, 1, 1), 0, 10)

