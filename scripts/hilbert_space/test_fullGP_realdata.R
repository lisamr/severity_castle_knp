# test out the hilbert space GP models on the fire CBI data

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
str(dat)
dat$preds %>% head
dat$preds$cbi %>% hist
summary(dat$preds$cbi)

# fix this later
dat$rings <- dat$rings %>% 
  mutate(dead_area.mean = dead_area.mean/900) # omit once you fix this


# assign an ordered, categorical variable for CBI
# cutpoints from Miller & Thode 2007 (unchanged, low, mod, high)
cuts <- c(-Inf, .1, 1.25, 2.25, Inf)
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



# remove NAs --------------------------------------------------------------

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_distance_dependence.R')

# all of the non-NA data
names(dat$preds)
non_NAs <- dat$preds %>% 
  select(cbi, ndvi, elevation:HLI, TPI_100:TPI_1000, count_red, count_grey) %>% 
  apply(2, function(x) !is.na(x)) %>% 
  apply(1, all) %>% which

# coords and response variables
x_pred_all <- dat$loc[non_NAs,c('x', 'y')]/1000
y_cat_all <- dat$preds$cbi_factor[non_NAs]

# density dependent variables (f_ddvar_wide zscores internally)
ddvar_wide_array <- 
  map(c('count_red.mean', 'count_grey.mean'), 
      ~f_ddvar_wide(dat$rings, .x)) %>% 
  map(~ .x[non_NAs,])

# other predictors
vars <- c('ndvi', 'HLI')
dat_preds <- dat$preds[non_NAs,] %>% 
  select(all_of(vars)) %>% 
  mutate(across(everything(), zscore))


# get area and distance vectors out
dist_vec <- dat$rings %>% 
  mutate(radius = radius / 1000) %>% 
  filter(id == 1) %>% 
  pull(radius) 
inner <- c(0, dist_vec[-length(dist_vec)]) # inner ring radius
area_vec <- f_ring_area(dist_vec, inner) 


# sample from larger dataset --------------------------------------------

# sample smaller set of data
set.seed(1)
N_sample <- 500
N_pred <- length(non_NAs)
idx <- sample(N_pred, N_sample)

# look at length scale prior
curve(pscl::densigamma(x, 2, .5), 0.01, 1) 
curve(pscl::densigamma(x, 4, .5), 0.01, 1, add = T, lty = 2) 
curve(dgamma(x, shape = 2, rate = 3), .01, 5)


# prep data list ----------------------------------------------------------

stan_fullGP_DD <- cmdstan_model('scripts/hilbert_space/stan_models/full_GP_ordinal_DD.stan')
print(stan_fullGP_DD)

dat_list <- list(
  N_sample = N_sample,
  N_pred = N_pred,
  idx = idx,
  nvars = ncol(dat_preds),
  K = max(y_cat_all),
  y = y_cat_all,
  coords = as.matrix(x_pred_all),
  X = as.matrix(dat_preds),
  nDD = length(ddvar_wide_array),
  R = length(dist_vec),
  dd_var = ddvar_wide_array,
  area = area_vec,
  dist_sq = dist_vec^2
)
str(dat_list)

# fit model
fitDD <- stan_fullGP_DD$sample(dat_list, iter_warmup = 1000, iter_sampling = 1000, 
                      refresh = 1, parallel_chains = 4,
                      output_dir = 'scripts/hilbert_space/model_outputs/',
                      output_basename = glue::glue("fullGP_DD_N{N_sample}"))



# other -------------------------------------------------------------------


# make data list
standata_full <- list(N_sample= N_sample,
                      N_pred= N_pred, 
                      vv_sample= idx, 
                      K = max(y_cat),
                      coords= x_pred, 
                      y_pred= y_cat)
str(standata_full)

stan_fullGP <- cmdstan_model('scripts/hilbert_space/stan_models/full_GP_ordinal.stan')

t1 <- proc.time()
fitGP <- stan_fullGP$sample(standata_full, parallel_chains = 4, refresh = 50, 
                            output_dir = 'scripts/hilbert_space/model_outputs/',
                            output_basename = glue::glue("fullGP_N{N_sample}"))
t2 <- proc.time()
t2 - t1
param <- c('rho', 'eta', 'cutpoint')
fitGP$summary(param)
# variable     mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
# <chr>       <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
#   1 rho          1.26   1.26 0.164 0.156  1.01   1.55   1.01     753.    1183.
# 2 eta          1.61   1.60 0.207 0.208  1.29   1.97   1.00    1331.    2488.
# 3 cutpoint[1] -4.32  -4.31 0.406 0.398 -5.00  -3.67   1.00    2781.    2849.
# 4 cutpoint[2] -1.20  -1.20 0.270 0.259 -1.65  -0.770  1.00    3000.    3420.
# 5 cutpoint[3]  1.07   1.07 0.268 0.267  0.641  1.52   1.00    3095.    3196.
mcmc_trace(fitGP$draws(c('rho', 'eta')))
dist(x_pred[idx,]) %>% c %>% hist(breaks = 100, xlim = c(0, 5))

# fit
y_rep <- fitGP$draws('y_rep', format = 'matrix')
y_sample <- with(standata_full, y_pred[vv_sample])
ppc_bars(y_sample, y_rep)




# look at gp
y_rep_mode <- apply(y_rep, 2, Mode)
GP_predictions <- x_pred[idx,] %>% 
  mutate(y_true = y_cat[idx],
         y_rep = y_rep_mode)
head(GP_predictions)

# y
cowplot::plot_grid(
  ggplot(GP_predictions, aes(x, y)) +
    geom_point(aes(color = y_true)) +
    coord_equal() +
    harrypotter::scale_color_hp(limits = c(0, 4)) +
    labs(title = 'true y'),
  ggplot(GP_predictions, aes(x, y)) +
    geom_point(aes(color = y_rep)) +
    coord_equal()+
    harrypotter::scale_color_hp(limits = c(0, 4)) +
    labs(title = 'posterior mode of y')
)
