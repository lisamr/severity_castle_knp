# run ordinal regressions on the CBI data with distance-dependent predictors
library(sf)
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
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# look at data
str(dat); names(dat)
dat$all_vars %>% head
dat$rings %>% head
dat$rings_conifer %>% head
dat$all_vars$cbi %>% hist
summary(dat$all_vars$cbi)
length(unique(dat$rings$id)) == nrow(dat$all_vars) # should be T



# small data manipulations ------------------------------------------------


# assign an ordered, categorical variable for CBI
# cutpoints from Miller & Thode 2007 (unchanged, low, mod, high)
cuts <- c(0, .1, 1.25, 2.25, 3)
dat$all_vars <- dat$all_vars %>% 
  mutate(cbi_factor = cut(cbi, cuts, labels = F, right = F, include.lowest = T))

# look at distribution of CBI
dat$all_vars %>% 
  ggplot(aes(cbi, fill = as.factor(cbi_factor))) +
  geom_histogram(binwidth = .01) +
  harrypotter::scale_fill_hp_d(option = 'Ravenclaw')

# merge rings together 
dat$rings <- full_join(st_drop_geometry(dat$rings), 
                       st_drop_geometry(dat$rings_conifer))

# get variables together --------------------------------------------------

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_distance_dependence.R')

f_prepdatlist <- function(X_vars, dd_vars, N, l1, l2,lscale_a, lscale_b, delta_a, delta_b, conifer_only = F, seed = NULL){
  # make sure dat$all_vars and $rings is in environment
  stopifnot(exists('dat'))
  stopifnot(all(c("all_vars", 'rings') %in% names(dat)))
  
  # prep data ===================
  
  # response variable
  y <- dat$all_vars$cbi_factor
  
  # controlling variables
  X <- dat$all_vars %>% 
    select(all_of(X_vars)) %>% 
    mutate(across(everything(), zscore)) %>% 
    as.matrix
  
  # need to check for and remove NAs
  stopifnot(sum(is.na(y)) == 0)
  stopifnot(sum(apply(X, 2, function(x) sum(is.na(x)))) == 0)
  
  # list of dist-dep variables
  ddvar_wide_array <- map(dd_vars, ~f_ddvar_wide(st_drop_geometry(dat$rings), .x))
  
  # get area and distance vectors out
  dist_vec <- dat$rings %>% 
    mutate(radius = radius / 1000) %>% 
    filter(id == 1) %>% 
    pull(radius) 
  inner <- c(0, dist_vec[-length(dist_vec)]) # inner ring radius
  area_vec <- f_ring_area(dist_vec, inner) 
  
  # put data into list ==========
  if(!is.null(seed)) {set.seed(seed)} 
  
  if(conifer_only){
    conifer_idx <- which(dat$all_vars$conifer == 1)
    idx <- sample(conifer_idx, N)
  }else{
    idx <- sample(length(y), N)
  }

  
  # deal with HSGP stuff ===========
  m_QE <- function(c,l,S) ceiling(1.75 * c / (l/S))
  coords <- dat$loc / 1000
  S1 <- (max(coords[,1]) - min(coords[,1]))/2
  S2 <- (max(coords[,2]) - min(coords[,2]))/2
  m1 <- m_QE(1.2, l1, S1)
  m2 <- m_QE(1.2, l2, S2)
  indices <- as.matrix(expand_grid(1:m1, 1:m2))
  
  # make a data list
  dat_list <- list(
    # HSGP stuff
    D = 2L,
    L = c(c(1.2, 1.2) * c(S1, S2)),
    M = c(m1, m2),
    M_nD = m1*m2,
    indices = indices,
    coords = coords,
    # regular stuff
    N = N,
    N_total = nrow(dat$loc),
    K = max(y),
    J = ncol(X),
    X = X,
    y = y,
    nDD = length(ddvar_wide_array),
    R = length(dist_vec),
    dd_var = ddvar_wide_array,
    area = area_vec,
    dist_sq = dist_vec,
    idx = idx,
    lscale_a = lscale_a,
    lscale_b = lscale_b,
    delta_a =delta_a,
    delta_b = delta_b
  )
  
  return(dat_list)
}




# prep data list ----------------------------------------------------------

# look at this for figuring out priors
curve(pscl::densigamma(x, 5, 10), 0.01, 5) 
curve(pscl::densigamma(x, 1, .1), 0.01, 1) 
curve(dgamma(x, shape = 2, rate = 5), .01, 5)
curve(dlnorm(x, -1, 1.5), .01, 2, ylim = c(0, 3))

# load up the model
HSGP <- cmdstan_model('scripts/stan_models/ordinal_HSGP.stan')
HSGP_DD <- cmdstan_model('scripts/stan_models/ordinal_HSGP_DD.stan')
HSGP_iso_DD <- cmdstan_model('scripts/stan_models/ordinal_isotropic_HSGP_DD.stan')

names(dat$all_vars)
names(dat$rings)
# "dead_area" "dead_count" "count_red" "count_grey" "TPA" "cover"         

vars <- c(
  'ndvi', 'TPA_conifer',
  # topography
  'elevation', 'slope', 'HLI', 'TPI_500',
  # climate
  'vpd_avg', #'wind_avg',
  # weather
  'vpd_day', 'wind_day'
)
#dd <- c('dead_area', 'cover')
dd <- c('dead_count_con')#, 'TPA_con') # same as area. 
dd <- c('dead_area_con', 'cover_con') # variables make sense. same results as counts
#dd <- c('count_red_con', 'count_grey_con', 'TPA_con')
N <- 1000

dat_list <- f_prepdatlist(vars, dd, N, l1 = 2, l2 = 2, 
                          lscale_a = 5, lscale_b = 10, delta_a = 1, delta_b = .1)
str(dat_list)

# fit the model
fit <- HSGP_DD$sample(dat_list, parallel_chains = 4, chains = 4, 
                   iter_sampling = 500, iter_warmup = 500, refresh = 1)
fit$summary(variables = c('B', 'B_dd', 'delta', 'cutpoint', 'lscale', 'gpscale'))
data.frame(id = 1:length(vars), vars); dd

mcmc_trace(fit$draws(c('cutpoint', 'lscale', 'gpscale', 'delta')))

# look at fit
y_rep <- fit$draws('y_rep', format = 'matrix')
ppc_bars(dat_list$y[dat_list$idx], y_rep)
#mcmc_pairs(fit$draws('B'))
#mcmc_pairs(fit$draws(c('B[1]', 'beta')))
#mcmc_pairs(fit$draws(c('beta', 'c', 'delta')))

# look at GP
gp <- fit$draws('f_pred', format = 'matrix')
p_gp <- as_tibble(dat$loc) %>% 
  mutate(gp_mean = colMeans(gp),
         gp1 = c(gp[1,])) %>% 
  ggplot(aes(x, y, fill = gp_mean)) +
  geom_raster() +
  coord_equal()

# look at predictions
y_rep_all <- fit$draws('y_rep_all', format = 'matrix')
y_rep_all_mode <- apply(y_rep_all, 2, Mode)
y_pred_df <- as_tibble(dat$loc) %>% 
  mutate(y_true = dat_list$y,
         y_rep_mode = y_rep_all_mode,
         y_rep1 = c(y_rep_all[1,])) %>% 
  mutate(across(y_true:y_rep1, as.factor))
p_yrep <- ggplot(y_pred_df, aes(x, y, fill = y_rep_mode)) +
  geom_raster() +
  coord_equal()+
  labs(title = 'posterior mode') +
  theme(legend.position = 'none')+
  harrypotter::scale_fill_hp_d('Ravenclaw')
p_yrep1 <- ggplot(y_pred_df, aes(x, y, fill = y_rep1)) +
  geom_raster() +
  coord_equal() +
  labs(title = 'first draw') +
  theme(legend.position = 'none')+
  harrypotter::scale_fill_hp_d('Ravenclaw')
p_ytrue <- ggplot(y_pred_df, aes(x, y, fill = y_true)) +
  geom_raster() +
  coord_equal()+
  labs(title = 'observed') +
  theme(legend.position = 'none')+
  harrypotter::scale_fill_hp_d('Ravenclaw')
cowplot::plot_grid(p_yrep, p_yrep1, p_ytrue, p_gp)


loglik <- fit$draws('log_lik', format = 'matrix')
loglik_mean <- colMeans(loglik)
as_tibble(dat$loc[dat_list$idx,]) %>% 
  mutate(loglik_mean) %>% 
  ggplot(aes(x, y, color = exp(loglik_mean))) +
  geom_point() +
  coord_equal() +
  scale_color_distiller(palette = 'RdBu')
