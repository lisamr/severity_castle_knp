# run ordinal regressions on the CBI data with distance-dependent predictors
library(sf)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(tidybayes)

rm(list = ls())

dat <- read_rds('outputs/tabular/model_knp_ringdata.rds')


theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
HPDI <- rethinking::HPDI
inv_logit <- rethinking::inv_logit
zscore <- function(x){
  stopifnot(is.numeric(x))
  # only do for non-dummy variables
  dummy <- all(unique(x[!is.na(x)]) %in% c(0,1))
  if(dummy) return(x) else{
    (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
  }
} 
zscore_sim <- function(xsim, x) (xsim - mean(x)) / sd(x) # for zscoring simulated data
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# for sampling on a smaller grid
idx_gridsmaller <- function(by){
  x_locs <- unique(dat$loc[,1]) %>% sort
  keep_x <- x_locs[seq(1, length(x_locs), by = by)]
  y_locs <- unique(dat$loc[,2]) %>% sort
  keep_y <- y_locs[seq(1, length(y_locs), by = by)]
  
  grid_sm <- dat$loc %>% 
    as_tibble() %>% 
    mutate(keep = case_when(
      x %in% keep_x & y %in% keep_y ~ 1, 
      T ~ 0) ) %>% 
    pull(keep) 
  
  which(grid_sm == 1)
}

# look at data
# str(dat); names(dat)
# dat$all_vars %>% head
# dat$rings %>% head
# dat$rings_conifer %>% head
# dat$all_vars$cbi %>% hist
# summary(dat$all_vars$cbi)
# length(unique(dat$rings$id)) == nrow(dat$all_vars) # should be T
# 


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
dat$rings2 <- full_join(st_drop_geometry(dat$rings), 
          st_drop_geometry(dat$rings_conifer)) #%>% filter(radius < 360)

# turn veg categories into a set of dummy variables
dat$all_vars <- dat$all_vars %>% 
  mutate(veg_2 = ifelse(veg == 2, 1, 0),
         veg_3 = ifelse(veg == 3, 1, 0),
         veg_4 = ifelse(veg == 4, 1, 0))



# get variables together --------------------------------------------------

source('scripts/functions/functions_process_spatial.R')
source('scripts/functions/functions_distance_dependence.R')

f_prepdatlist <- function(X_vars, dd_vars, int_var, .N=NULL, .by = NULL, 
                          mask_idx = NULL, seed = NULL){
  # make sure dat$all_vars and $rings is in environment
  stopifnot(exists('dat'))
  stopifnot(all(c("all_vars", 'rings2') %in% names(dat)))
  if(is.null(.N) & is.null(.by)) stop('either N or by needs to be specified.')
  if(!is.null(.N) & !is.null(.by)) warning("N will by overridden by 'by'")
  
  
  # select data =================
  # response variable
  y <- dat$all_vars$cbi_factor
  
  # controlling variables
  X <- dat$all_vars %>% 
    select(all_of(X_vars))
  
  # list of dist-dep variables
  ddvar_wide_array <- map(dd_vars, ~f_ddvar_wide(st_drop_geometry(dat$rings2), .x))
  
  # identify masked rows ==============
  # NAs
  NA_y <- which(is.na(y))
  NA_X <- unname(which(apply(X, 1, function(x) any(is.na(x)))))
  NAs <- c(NA_y, NA_X)
  
  # anything else
  mask_idx <- c(NAs, mask_idx)
  

  # sample data to smaller subset =====
  if(!is.null(seed)) {set.seed(seed)} 
  
  # if you want to sample on a regular grid
  if(!is.null(.by)){
    idx <- idx_gridsmaller(.by)
    # remove any masked points
    idx <- idx[!idx %in% mask_idx]
    .N <- length(idx)
  }else{ 
    # if just sample N
    stopifnot(!is.null(.N) & is.null(.by))
    idx <- seq_along(y)[!seq_along(y) %in% mask_idx]
    idx <- sample(idx, .N)
  }
  
  
  # prep data =========================
  
  # zscore controlling variables
  X <- apply(X, 2, zscore)

  # get area and distance vectors out
  dist_vec <- dat$rings2 %>% 
    mutate(radius = radius / 1000) %>% 
    filter(id == 1) %>% 
    pull(radius) 
  inner <- c(0, dist_vec[-length(dist_vec)]) # inner ring radius
  area_vec <- f_ring_area(dist_vec, inner) 
  
  
  
  
  # put data into list ==========
  

  # make a data list
  dat_list <- list(
    N = .N,
    K = max(y[idx]),
    J = ncol(X),
    X = X[idx,],
    y = y[idx],
    nDD = length(ddvar_wide_array),
    M = length(dist_vec),
    dd_var = map(ddvar_wide_array, ~ .x[idx,]),
    area = area_vec,
    dist_sq = dist_vec^2,
    idx = idx,
    nInt = length(int_var),
    int_var = match(int_var, colnames(X))
  )
  
  return(dat_list)
}

# make function for adding simulation data into data list

# create data to simulate from 
f_add_simdata <- function(dat_list, ddvar, intvar){
  
  stopifnot(length(ddvar) == length(intvar))
  
  dd_var_sim <- matrix(rep(ddvar, dat_list$M), ncol = dat_list$M)
  Xsim <- intvar
  Nsim <- nrow(dd_var_sim)
  # add to dat_list
  dat_list$dd_var_sim <- dd_var_sim
  dat_list$Xsim <- Xsim
  dat_list$Nsim <- Nsim
  return(dat_list)
}




# prep data list ----------------------------------------------------------

# load up the model
stanmod <- cmdstan_model('scripts/stan_models/ordinal_DDarray.stan')
stanmod_int <- cmdstan_model('scripts/stan_models/ordinal_DDarray_int.stan')

names(dat$all_vars)
names(dat$rings2)
# "dead_area" "dead_count" "count_red" "count_grey" "TPA" "cover"         

vars <- c(
  # vegetation
  'ndvi', 
  'veg_2', 'veg_3', #'TPA_conifer', # intercept=conifer, # now: 1 = conifer, veg2 = hardwood, veg3 = shrub/herbaceous
  # topography
  'elevation', 'slope', 'HLI', 'TPI_500',
  # climate
  'vpd_avg', #'wind_avg',
  # weather
  'vpd_day'#, 'wind_day'
)
#dd <- c('dead_area', 'cover') # tree density without masking out non-conifer vegetation
#dd <- c('dead_area_con', 'cover_con') # tree density (area) with non-conifers masked out.
dd <- c('dead_count_con', 'TPA_con') # tree density (counts) with non-conifers masked out
intvars <- c('vpd_day', 'veg_2', 'veg_3')#, 'wind_day')
N <- 5000
by <- 2 # every n cells sampled on a regular grid
mask_ids <- which(dat$all_vars$veg == 4) # mask out veg4

dat_list <- f_prepdatlist(X_vars = vars, dd_vars = dd, int_var = intvars,
                          .N = N, 
                          .by = by, 
                          mask_idx = mask_ids)
# add in data to simulate
#dat$all_vars$dead_count_con %>% hist(breaks = 15)
vpd_quants <- quantile(dat$all_vars$vpd_day, probs = c(.1, .5, .9))
sim_grid <- expand_grid(ddvar = seq(0,5,by=.2),
                         intvar = vpd_quants) %>% 
  mutate(ddvarz =  zscore_sim(ddvar, dat$rings2$dead_count_con),
         intvarz = unname(zscore_sim(intvar, dat$all_vars$vpd_day))
         )
dat_list <- f_add_simdata(dat_list, sim_grid$ddvarz, sim_grid$intvarz)
str(dat_list)

# see what you've sampled. 
tmp <- as.data.frame(dat$loc[dat_list$idx,])
ggplot(tmp, aes(x, y)) +
  geom_point() +
  coord_equal()
spacing <- min(dist(tmp[1:10,]))
spacing # spacing of grid in meters


# fit the model
#fit <- stanmod$sample(dat_list, parallel_chains = 4, chains = 4, iter_sampling = 1000, iter_warmup = 1000)
#fit$summary(variables = c('B', 'beta', 'c', 'delta'))
fit_int <- stanmod_int$sample(dat_list, parallel_chains = 4, chains = 4, iter_sampling = 1000, iter_warmup = 1000)
fit_int$summary(variables = c('B', 'beta', 'beta_int', 'c', 'delta')) %>% print(n = Inf)
data.frame(id = 1:length(vars), vars); dd


# look at fit
y_rep <- fit_int$draws('y_rep', format = 'matrix')
ppc_bars(dat_list$y, y_rep)
#mcmc_pairs(fit$draws('B'))
#mcmc_pairs(fit$draws(c('B[1]', 'beta')))
#mcmc_pairs(fit$draws(c('beta', 'c', 'delta')))

# log likelihood
loglik <- fit$draws('log_lik', format = 'matrix')
loglik_mean <- colMeans(loglik)
yrep_mode <- c(apply(y_rep, 2, Mode))

errors_df <- tibble(
  x = dat$loc[dat_list$idx,1],
  y = dat$loc[dat_list$idx,2],
  #loglik_mean, 
  cbi = dat_list$y,
  yrep_mode,
  y_rep1 = c(y_rep[1,]))
p_loglik <- ggplot(errors_df, aes(x, y, color = (loglik_mean))) +
  geom_point() +
  coord_equal()+
  scale_color_distiller(palette = 'RdBu')
p_y <- ggplot(errors_df, aes(x, y, color = as.factor(cbi))) +
  geom_point() +
  coord_equal() +
  scale_color_brewer(palette = 'RdBu', direction = -1) + 
  theme(legend.position = 'bottom')
p_yrep <- ggplot(errors_df, aes(x, y, color = (y_rep1))) +
  geom_point() +
  coord_equal()+
  scale_color_distiller(palette = 'RdBu')+ 
  theme(legend.position = 'bottom')
p_yrepm <- ggplot(errors_df, aes(x, y, color = as.factor(yrep_mode))) +
  geom_point() +
  coord_equal()+
  scale_color_brewer(palette = 'RdBu', direction = -1, limits = factor(1:4))+ 
  theme(legend.position = 'bottom')

cowplot::plot_grid(
  p_loglik, 
  p_y,
  p_yrep,
  p_yrepm
)
cowplot::plot_grid(
  p_y,
  p_yrepm
)

# look at correlations between variables
dat_list
dat$all_vars[dat_list$idx,c('elevation', 'vpd_avg', 'rh_avg', 'wind_avg')] %>% pairs
cor(dat$all_vars[dat_list$idx,c('elevation', 'vpd_avg', 'rh_avg', 'wind_avg')])
dat$all_vars[dat_list$idx,c('elevation', 'vpd_day', 'rh_day', 'wind_day')] %>% pairs
cor(dat$all_vars[dat_list$idx,c('elevation', 'vpd_day', 'rh_day', 'wind_day')])


dat$all_vars[dat_list$idx,c('dead_count_con', 'TPA_conifer', 'elevation', 'vpd_avg', 'vpd_day')] %>% pairs



# check out simulated data ------------------------------------------------

sim_grid
phi_sim <- fit_int$draws(c('phi_sim'), format = 'matrix')
draws_c <- fit_int$draws(c('c'), format = 'matrix')

# cumulative proportions, in list form. 
x <- lapply(1:ncol(draws_c), function(k){
  # equivalent to: inv_logit(c - phi)
  inv_logit(sweep(-phi_sim, 1, draws_c[,k], FUN = '+')) 
})

# get posterior proportions for simulated data. 
post_sim <- lapply(seq_len(ncol(phi_sim)), 
                   # combine the columns of each matrix across the list
                   function(i) do.call(cbind, lapply(x, '[', , i))) %>% 
  # convert to proportions in each CBI class
  map(~ cbind(0, .x, 1)) %>% 
  map(~ apply(.x, 1, diff) %>% t) %>% 
  map(~ tibble(class = 1:(ncol(.x)),
               median = apply(.x, 2, median),
               lower = apply(.x, 2, HPDI, .9)[1,],
               upper = apply(.x, 2, HPDI, .9)[2,])) %>% 
  map2(.x = ., .y = 1:length(.), ~ cbind(.x, index = .y)) %>% 
  bind_rows() %>% 
  # merge with dat_sim
  left_join(sim_grid %>% mutate(index = 1:nrow(.))) %>% 
  mutate(class = as.factor(class),
         dead_TPH = ddvar/900*10000)

# class palette
pal <- wesanderson::wes_palette('Zissou1', n = 15, type = 'continuous')
pal
Pal <- pal[c(1, 7, 12, 15)]
labelnames <- c('10% VPD', '50% VPD', '90% VPD')
names(labelnames) <- unique(sim_grid$intvar)
classlabels <- c('unchanged', 'low', 'moderate', 'high')
names(classlabels) <- 1:4

ggplot(post_sim, aes(dead_TPH, median, group = class, color = class, fill = class)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha= .6, color = NA) +
  facet_wrap(~intvar, labeller = labeller(intvar = labelnames)) +
  labs(x = 'Dead trees/ha', y = 'Probability of severity class', 
       color = 'Severity class', fill = 'Severity class') +
  scale_color_manual(values = Pal, labels = classlabels) +
  scale_fill_manual(values = Pal, labels = classlabels) +
  theme(legend.position = 'bottom') 
ggsave('figures/deadcounts_vpd.png', width = 8, height = 4)

dat$all_vars %>% 
  mutate(dead_TPH = dead_count_con/900*10000,
         cbi_factor = as.factor(cbi_factor)) %>% 
  ggplot(aes(dead_TPH, group = cbi_factor, fill = cbi_factor)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 50)) +
  facet_grid(~ cbi_factor, labeller = labeller(cbi_factor = classlabels)) +
  labs(x = 'Dead trees/ha', y = 'density', 
       color = 'Severity class', fill = 'Severity class') +
  scale_fill_manual(values = Pal, labels = classlabels) +
  theme(legend.position = 'bottom') 
ggsave('figures/deadcounts.png', width = 8, height = 4)





# posterior estimates of kernels ------------------------------------------


# generate distance and area vectors to simulate w to new data

delta <- fit_int$draws('delta[1]', format = 'matrix')
dist_vec_sim <- seq(.05, 1, by = .01)
f_simulatew <- function(delta, dist_vec_sim){
  dist_vec <- sqrt(dat_list$dist_sq)
  area_vec <- dat_list$area
  lm_a_d <- lm(area_vec ~ dist_vec) # luckily there's a linear relationship
  area_vec_sim <- predict(lm_a_d, list(dist_vec = dist_vec_sim))
  
  sapply(1:nrow(delta), function(i){
    w0 <- exp(-dist_vec_sim^2 / (2*(delta[i]^2)) ) * area_vec_sim
    w <- w0/sum(w0)
    return(w)
  }) %>% t
} 

w <- f_simulatew(delta, dist_vec_sim)


#png('figures/kernels.png', width = 8, height = 4, units = 'in', res = 300)
par(mfrow = c(1,2))
f_plot_cov(delta, seq(0, 1, by = .01), xlim = c(.05, .4), 
           main = 'Decay kernel of dead trees', xlab = 'Distance (km)', ylab = 'Correlation')
f_w90(w, xlim= c(.05, .4), xlab = 'Distance (km)', ylab = 'Cumulative effect', 
      main = 'Cumulative effect of dead trees')
#dev.off()




# plotting covariates -----------------------------------------------------

dd_names <- c('dead.trees', 'live.trees')
fit_summary <- fit_int$summary(variables = c('B', 'beta', 'c', 'delta')) %>% 
  mutate(variable_name = c(vars, dd_names, 
                           #paste0('vpd x ', dd_names), 
                           paste0("c_", 1:3), 
                           paste0("delta_", dd_names)), 
         variable_namef = factor(variable_name, levels = rev(variable_name)),
         dashed = ifelse(q5*q95 > 0, F, T),
         color = ifelse(mean > 0 , "#F21A00", "#3B9AB2"),
         .after = variable) 
print(fit_summary, n = Inf)

fit_summary2 <- fit_summary %>% 
  filter(!str_detect(variable_namef, 'c_'), 
         !str_detect(variable_namef, 'delta_'))
  
ggplot(fit_summary2, aes(mean, variable_namef)) +
  geom_vline(xintercept = 0, linewidth = 4, color = 'grey90') +
  geom_pointrange(aes(xmin = q5, xmax = q95, linetype = dashed), 
                  color = fit_summary2$color, linewidth = 1, size = 1) +
  labs(x = 'Log-odds (90% HPDI)', y = 'Parameters') +
  theme(legend.position = 'none')

ggsave('figures/covariates.png', width = 6, height = 4)
