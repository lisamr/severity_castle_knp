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



# sample from larger dataset --------------------------------------------

# all of the non-NA data
non_NAs <- which(!is.na(dat$preds$cbi)) 
x_pred_all <- dat$loc[non_NAs,c('x', 'y')]/1000
y_cat_all <- dat$preds$cbi_factor[non_NAs]

# tune the model on a much smaller part of the dataset
xrange <- c(335, 340)
yrange <- c(4045, 4055)
subset_data <- x_pred_all %>% 
  mutate(y_cat = y_cat_all) %>% 
  filter(x >= xrange[1], x <= xrange[2],
         y >= yrange[1], y <= yrange[2])
hist(subset_data$y_cat)

ggplot(x_pred_all %>% mutate(y_cat = y_cat_all), aes(x, y)) +
  geom_point() +
  geom_point(data = subset_data, aes(color = y_cat), size = 2) +
  coord_equal() +
  harrypotter::scale_color_hp(option = 'HarryPotter', direction = -1)

# sample data for tuning
set.seed(1)
N_sample <- 500

# prep data for model
x_pred <- select(subset_data, x, y)
y_cat <- subset_data$y_cat

# sample non-NA data
N_pred <- length(y_cat)
idx <- sample(N_pred, N_sample)


# Half range S for each dimension
S1 <- abs(max(x_pred$x) - min(x_pred$x)) / 2
S2 <- abs(max(x_pred$y) - min(x_pred$y)) / 2
Dims <- ncol(x_pred)

# look at length scale prior
curve(pscl::densigamma(x, 2, .5), 0.01, 1) 
curve(pscl::densigamma(x, 4, .5), 0.01, 1, add = T, lty = 2) 
curve(dgamma(x, shape = 2, rate = 3), .01, 5)



# relationships among m, l, c ---------------------------------------------

m_QE <- function(c,l,S) ceiling(1.75 * c / (l/S))
l_QE <- function(c,m,S) round(S * 1.75 * c / m, 3)   
c_vs_l_QE <- function(l,S){
  c =  3.2*l/S
  if(c < 1.2)
    c = 1.2
  c
}
diagnostic <- function(l,l_hat) l_hat + 0.01 >= l

standata <- list()
fit<- list()
diagnosis <- list()

l1 <- vector()  # lengthscale for dimension 1
l2 <- vector()  # lengthscale for dimension 2
c1 <- vector()  # boundary factor for dimension 1
c2 <- vector()  # boundary factor for dimension 2
m1 <- vector()  # number of basis functions for dimension 1
m2 <- vector()  # number of basis functions for dimension 2
l_hat1 <- vector()  # HSGP lengthscale estimate for dimension 1
l_hat2 <- vector()  # HSGP lengthscale estimate for dimension 2
check1 <- vector()  # Diagnostic for dimension 1
check2 <- vector()  # Diagnostic for dimension 2
rmse_gp <- vector()  # root mean square error of GP model
rmse <- vector()  # root mean square error of HSGP model
eR <- vector()  # coefficient of determination of HSGP model
elpd <- vector()  # expected log predictive density of HSGP model


stan_ordinal1 <- cmdstan_model('scripts/hilbert_space/stan_models/HSGP_2D_ordinal_novars.stan')
stan_isotropic1 <- cmdstan_model('scripts/hilbert_space/stan_models/HSGP_2D_istotropic.stan')


# tune parameters: isotropic -------------------------------------------------

# setting m, l, c
i = 1
runmod <- F
imax <- 1
for(i in 1:imax){
  print(paste("Iteration i =", i), quote = FALSE)
  if(i==1){
    l1[i] <- 1
    c1[i] <- c_vs_l_QE(l=l1[i], S=S1)
    m1[i] <- m_QE(c=c1[i], l=l1[i], S=S1)
    
    l2[i] <- 1
    c2[i] <- c_vs_l_QE(l=l2[i], S=S2)
    m2[i] <- m_QE(c=c2[i], l=l2[i], S=S2)
  }else{
    # dim1
    if(check1[i-1]){
      if(i > 2){
        if(check1[i-2]){ # T, T
          m1[i] <- m1[i-1]
          c1[i] <- c1[i-1]
          l1[i] <- l1[i-1]
        }else{ # T, F
          m1[i] <- m1[i-1] + 2
          c1[i] <- c_vs_l_QE(l=l_hat1[i-1], S=S1)
          l1[i] <- l_QE(c=c1[i], m=m1[i], S=S1)
        }
      }else{
        m1[i] <- m1[i-1] + 2
        c1[i] <- c_vs_l_QE(l=l_hat1[i-1], S=S1)
        l1[i] <- l_QE(c=c1[i], m=m1[i], S=S1)
      }
    }else{ # F
      l1[i] <- l_hat1[i-1]
      c1[i] <- c_vs_l_QE(l=l1[i], S=S1)
      m1[i] <- m_QE(c=c1[i], l=l1[i], S=S1)
    }
    # dim2
    if(check2[i-1]){
      if(i > 2){
        if(check2[i-2]){ # T, T
          m2[i] <- m2[i-1]
          c2[i] <- c2[i-1]
          l2[i] <- l2[i-1]
        }else{ # T, F
          m2[i] <- m2[i-1] + 2
          c2[i] <- c_vs_l_QE(l=l_hat2[i-1], S=S2)
          l2[i] <- l_QE(c=c2[i], m=m2[i], S=S2)
        }
      }else{
        m2[i] <- m2[i-1] + 2
        c2[i] <- c_vs_l_QE(l=l_hat2[i-1], S=S2)
        l2[i] <- l_QE(c=c2[i], m=m2[i], S=S2)
      }
    }else{ # F
      l2[i] <- l_hat2[i-1]
      c2[i] <- c_vs_l_QE(l=l2[i], S=S2)
      m2[i] <- m_QE(c=c2[i], l=l2[i], S=S2)
    }
  }
  print(c(l1 = l1, c1 = c1, m1 = m1, S1 = S1) )
  print(c(l2 = l2, c2 = c2, m2 = m2, S2 = S2))
  
  
  # D-tuples
  indices <- as.matrix(expand_grid(1:m1[i], 1:m2[i]))
  
  # make data list
  standata[[i]] <- list(D= Dims, # dimensions
                        M= c(m1[i], m2[i]), # no. basis functions
                        M_nD= m1[i]*m2[i], # length of tuples
                        indices= indices, # D-tuples
                        L= c(c1[i], c2[i]) * c(S1, S2), # boundaries
                        x_pred= x_pred, # predictors for whole dataset
                        K = max(y_cat),
                        y_pred= y_cat, # response for whole dataset
                        N_sample= N_sample, # N sample
                        N_pred= N_pred, # N for whole dataset
                        vv_sample= idx # sampled indices
  )
  str(standata[[i]])
  
  if(runmod){
    # fit model
    t1 <- proc.time()
    fit[[i]] <- stan_isotropic1$sample(data= standata[[i]],
                                     iter_warmup=100,
                                     iter_sampling=100, refresh = 20,
                                     chains=4, parallel_chains = 4, 
                                     thin=2, init=0.5,
                                     # adapt_delta=0.9,
                                     save_warmup=TRUE)
    t2 <- proc.time()
    tdelta <- t2 - t1
    print(tdelta)
    
    # model evaluation
    f_hsgp <- fit[[i]]$summary('f')
    
    # rmse against data
    residual <- standata[[i]]$y_pred[idx] - f_hsgp$mean
    rmse[i] <- round(sqrt(mean(residual^2)), 2)
    
    # elpd
    elpd[i] <- round(median(fit[[i]]$summary('elpd')$mean), 2)
    
    l_hat1[i] <- fit[[i]]$summary('lscale')$mean[1]
    l_hat2[i] <- fit[[i]]$summary('lscale')$mean[1] #fit[[i]]$summary('lscale')$mean[2]
    check1[i] <- diagnostic(l1[i], l_hat1[i])
    check2[i] <- diagnostic(l2[i], l_hat2[i])
    
    
    diagnosis[[i]] <- data.frame(iter= c(i, i),
                                 Dim= c("D1", "D2"),
                                 l= c(l1[i], l2[i]),
                                 c= c(c1[i], c2[i]),
                                 m= c(m1[i], m2[i]),
                                 l_hat= c(l_hat1[i], l_hat2[i]),
                                 'l_hat > l' = c(check1[i], check2[i]),
                                 rsme = c(NA, rmse[i]),
                                 elpd = c(NA, elpd[i]), 
                                 tdelta = unname(tdelta[3]))
    names(diagnosis[[i]]) <- c("iter", "Dim", "l", "c", "m", "l_hat", "l_hat + 0.01 > l", 
                               "rmse", "elpd", 'tdelta') 
    
    if(i==1){ 
      print(diagnosis[[i]])
    }else{
      diagnosis[[i]] <- rbind(diagnosis[[i-1]],diagnosis[[i]])
      print(diagnosis[[i]])
    }
  }
  
}


# do for larger set of data -----------------------------------------------

# sample data for tuning
set.seed(1)
N_sample <- 2000

# prep data for model
x_pred <- x_pred_all
y_cat <- y_cat_all

# sample non-NA data
N_pred <- length(y_cat)
idx <- sample(N_pred, N_sample)


# Half range S for each dimension
S1 <- abs(max(x_pred$x) - min(x_pred$x)) / 2
S2 <- abs(max(x_pred$y) - min(x_pred$y)) / 2
Dims <- ncol(x_pred)


# iter Dim         l        c  m     l_hat l_hat + 0.01 > l rmse  elpd tdelta
# 1    1  D1 1.0000000 1.285141  6 0.8624282            FALSE   NA    NA   9.08
# 2    1  D2 1.0000000 1.200000 11 0.8624282            FALSE 2.87 -0.61   9.08
# 3    2  D1 0.8624282 1.200000  7 0.5473878            FALSE   NA    NA   9.59
# 4    2  D2 0.8624282 1.200000 13 0.5473878            FALSE 3.06 -0.76   9.59
# 5    3  D1 0.5473878 1.200000 10 0.6146146             TRUE   NA    NA  14.19
# 6    3  D2 0.5473878 1.200000 20 0.6146146             TRUE 3.06 -0.74  14.19
# 7    4  D1 0.4360000 1.200000 12 0.5715072             TRUE   NA    NA  26.83
# 8    4  D2 0.4750000 1.200000 22 0.5715072             TRUE 3.05 -0.74  26.83

{
  l1_full <- .57#diagnosis[[4]]$l_hat %>% tail(1)
  c1_full <- c_vs_l_QE(l=l1_full, S=S1)
  m1_full <- m_QE(c=c1_full, l=l1_full, S=S1)
  
  l2_full <- l1_full
  c2_full <- c_vs_l_QE(l=l2_full, S=S2)
  m2_full <- m_QE(c=c2_full, l=l2_full, S=S2)
  
  print(c(l1 = l1_full, c1 = c1_full, m1 = m1_full, S1 = S1) )
  print(c(l2 = l2_full, c2 = c2_full, m2 = m2_full, S2 = S2))
  
  
  # D-tuples
  indices <- as.matrix(expand_grid(1:m1_full, 1:m2_full))
  
  # make data list
  standata_full <- list(D= Dims, # dimensions
                        M= c(m1_full, m2_full), # no. basis functions
                        M_nD= m1_full*m2_full, # length of tuples
                        indices= indices, # D-tuples
                        L= c(c1_full, c2_full) * c(S1, S2), # boundaries
                        x_pred= x_pred, # predictors for whole dataset
                        K = max(y_cat),
                        y_pred= y_cat, # response for whole dataset
                        N_sample= N_sample, # N sample
                        N_pred= N_pred, # N for whole dataset
                        vv_sample= idx # sampled indices
  )
}
str(standata_full)



# check out final model ---------------------------------------------------



#fit_final <- fit[[i]]
t1 <- proc.time()
fit_final <- stan_isotropic1$sample(data= standata_full,
                                    refresh = 1,
                                  iter_warmup=1000,
                                  iter_sampling=1000,
                                  chains=4, init=0.5, parallel_chains = 4,
                                  output_dir = 'scripts/hilbert_space/model_outputs/',
                                  output_basename = glue::glue("HSGP_N{N_sample}")
)
t2 <- proc.time()
tdelta <- t2 - t1
print(tdelta)
# user  system elapsed 
# 122.45    7.94 7965.78 

#6:50pm-7:30. 40 mins.
timestamp()

param <- c('lscale', 'gpscale', 'cutpoint')
fit_final$summary(param)
# variable      mean median    sd   mad     q5    q95  rhat ess_bulk ess_tail
# <chr>        <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>    <dbl>    <dbl>
# 1 lscale       1.38   1.36  0.258 0.234  1.02   1.86   1.00     969.    1405.
# 2 gpscale      1.62   1.59  0.277 0.261  1.22   2.11   1.00    1776.    2476.
# 3 cutpoint[1] -4.17  -4.16  0.311 0.314 -4.68  -3.68   1.00    3321.    2698.
# 4 cutpoint[2] -1.22  -1.21  0.139 0.138 -1.45  -0.992  1.00    4661.    3349.
# 5 cutpoint[3]  0.838  0.837 0.134 0.134  0.621  1.06   1.00    3588.    3645.
mcmc_trace(fit_final$draws(c('lscale', 'gpscale')))

# fit
y_rep <- fit_final$draws('y_rep', format = 'matrix')
y_sample <- with(standata_full, y_pred[vv_sample])
ppc_bars(y_sample, y_rep)

# look at gp
f_pred <- fit_final$draws('f_pred', format = 'matrix')
y_rep_all <- fit_final$draws('y_rep_all', format = 'matrix')
f_pred_mean <- colMeans(f_pred)
y_rep_mode <- apply(y_rep_all, 2, Mode)

n <- sample(4000, 1)
GP_predictions <- x_pred %>% 
  mutate(f_pred = f_pred_mean,
         y_true = y_cat,
         y_rep = y_rep_mode,
         y_rep_n = as.numeric(y_rep_all[n,]))


# y
cowplot::plot_grid(
  ggplot(GP_predictions, aes(x, y)) +
    geom_point(aes(color = y_true), size = 1) +
    coord_equal() +
    harrypotter::scale_color_hp(limits = c(0, 4)) +
    labs(title = 'true y'),
  ggplot(GP_predictions, aes(x, y)) +
    geom_point(aes(color = y_rep_n), size = 1) +
    coord_equal()+
    harrypotter::scale_color_hp(limits = c(0, 4)) +
    labs(title = 'posterior mode of y')
)

ggsave('ypredicted_ordinal.png', width = 6, height = 4)

# GP
ggplot(GP_predictions, aes(x, y)) +
  geom_point(aes(color = f_pred)) +
  coord_equal() +
  harrypotter::scale_color_hp(limits = c(-4, 4)) +
  labs(title = 'estimated GP')



# tune parameters ---------------------------------------------------------

# setting m, l, c
i = 1
runmod <- F
imax <- 1
for(i in 1:imax){
  print(paste("Iteration i =", i), quote = FALSE)
  if(i==1){
    l1[i] <- 1
    c1[i] <- c_vs_l_QE(l=l1[i], S=S1)
    m1[i] <- m_QE(c=c1[i], l=l1[i], S=S1)
    
    l2[i] <- 1
    c2[i] <- c_vs_l_QE(l=l2[i], S=S2)
    m2[i] <- m_QE(c=c2[i], l=l2[i], S=S2)
  }else{
    # dim1
    if(check1[i-1]){
      if(i > 2){
        if(check1[i-2]){ # T, T
          m1[i] <- m1[i-1]
          c1[i] <- c1[i-1]
          l1[i] <- l1[i-1]
        }else{ # T, F
          m1[i] <- m1[i-1] + 2
          c1[i] <- c_vs_l_QE(l=l_hat1[i-1], S=S1)
          l1[i] <- l_QE(c=c1[i], m=m1[i], S=S1)
        }
      }else{
        m1[i] <- m1[i-1] + 2
        c1[i] <- c_vs_l_QE(l=l_hat1[i-1], S=S1)
        l1[i] <- l_QE(c=c1[i], m=m1[i], S=S1)
      }
    }else{ # F
      l1[i] <- l_hat1[i-1]
      c1[i] <- c_vs_l_QE(l=l1[i], S=S1)
      m1[i] <- m_QE(c=c1[i], l=l1[i], S=S1)
    }
    # dim2
    if(check2[i-1]){
      if(i > 2){
        if(check2[i-2]){ # T, T
          m2[i] <- m2[i-1]
          c2[i] <- c2[i-1]
          l2[i] <- l2[i-1]
        }else{ # T, F
          m2[i] <- m2[i-1] + 2
          c2[i] <- c_vs_l_QE(l=l_hat2[i-1], S=S2)
          l2[i] <- l_QE(c=c2[i], m=m2[i], S=S2)
        }
      }else{
        m2[i] <- m2[i-1] + 2
        c2[i] <- c_vs_l_QE(l=l_hat2[i-1], S=S2)
        l2[i] <- l_QE(c=c2[i], m=m2[i], S=S2)
      }
    }else{ # F
      l2[i] <- l_hat2[i-1]
      c2[i] <- c_vs_l_QE(l=l2[i], S=S2)
      m2[i] <- m_QE(c=c2[i], l=l2[i], S=S2)
    }
  }
  print(c(l1 = l1, c1 = c1, m1 = m1, S1 = S1) )
  print(c(l2 = l2, c2 = c2, m2 = m2, S2 = S2))
  
  
  # D-tuples
  indices <- as.matrix(expand_grid(1:m1[i], 1:m2[i]))
  
  # make data list
  standata[[i]] <- list(D= Dims, # dimensions
                        M= c(m1[i], m2[i]), # no. basis functions
                        M_nD= m1[i]*m2[i], # length of tuples
                        indices= indices, # D-tuples
                        L= c(c1[i], c2[i]) * c(S1, S2), # boundaries
                        x_pred= x_pred, # predictors for whole dataset
                        K = max(y_cat),
                        y_pred= y_cat, # response for whole dataset
                        N_sample= N_sample, # N sample
                        N_pred= N_pred, # N for whole dataset
                        vv_sample= idx # sampled indices
  )
  str(standata[[i]])
  
  if(runmod){
    # fit model
    t1 <- proc.time()
    fit[[i]] <- stan_ordinal1$sample(data= standata[[i]],
                                    iter_warmup=100,
                                    iter_sampling=100, refresh = 10,
                                    chains=4, parallel_chains = 4, 
                                    thin=2, init=0.5,
                                    # adapt_delta=0.9,
                                    save_warmup=TRUE)
    t2 <- proc.time()
    tdelta <- t2 - t1
    print(tdelta)
    
    # model evaluation
    f_hsgp <- fit[[i]]$summary('f')
    
    # rmse against data
    residual <- standata[[i]]$y_pred[idx] - f_hsgp$mean
    rmse[i] <- round(sqrt(mean(residual^2)), 2)
    
    # elpd
    elpd[i] <- round(median(fit[[i]]$summary('elpd')$mean), 2)
    
    l_hat1[i] <- fit[[i]]$summary('lscale')$mean[1]
    l_hat2[i] <- fit[[i]]$summary('lscale')$mean[2]
    check1[i] <- diagnostic(l1[i], l_hat1[i])
    check2[i] <- diagnostic(l2[i], l_hat2[i])
    
    
    diagnosis[[i]] <- data.frame(iter= c(i, i),
                                 Dim= c("D1", "D2"),
                                 l= c(l1[i], l2[i]),
                                 c= c(c1[i], c2[i]),
                                 m= c(m1[i], m2[i]),
                                 l_hat= c(l_hat1[i], l_hat2[i]),
                                 'l_hat > l' = c(check1[i], check2[i]),
                                 rsme = c(NA, rmse[i]),
                                 elpd = c(NA, elpd[i]), 
                                 tdelta = unname(tdelta[3]))
    names(diagnosis[[i]]) <- c("iter", "Dim", "l", "c", "m", "l_hat", "l_hat + 0.01 > l", 
                               "rmse", "elpd", 'tdelta') 
    
    if(i==1){ 
      print(diagnosis[[i]])
    }else{
      diagnosis[[i]] <- rbind(diagnosis[[i-1]],diagnosis[[i]])
      print(diagnosis[[i]])
    }
  }
  
}

# 5:21
timestamp()


# check out final model ---------------------------------------------------

str(standata[[i]])

#fit_final <- fit[[i]]

fit_final <- stan_ordinal1$sample(data= standata[[i]],
                                 iter_warmup=1000,
                                 iter_sampling=1000,
                                 chains=4, init=0.5, parallel_chains = 4
)
#12:18
param <- c('lscale', 'gpscale', 'cutpoint')
fit_final$summary(param)
mcmc_trace(fit_final$draws(c('lscale', 'gpscale')))

# fit
y_rep <- fit_final$draws('y_rep', format = 'matrix')
y_sample <- with(standata[[i]], y_pred[vv_sample])
ppc_bars(y_sample, y_rep)

# look at gp
f_pred <- fit_final$draws('f_pred', format = 'matrix')
y_rep_all <- fit_final$draws('y_rep_all', format = 'matrix')
f_pred_mean <- colMeans(f_pred)
y_rep_mode <- apply(y_rep_all, 2, Mode)
GP_predictions <- x_pred %>% 
  mutate(f_pred = f_pred_mean,
         y_true = y_cat,
         y_rep = y_rep_mode)
head(GP_predictions)

# GP
ggplot(GP_predictions, aes(x, y)) +
  geom_point(aes(color = f_pred)) +
  coord_equal() +
  harrypotter::scale_color_hp(limits = c(-4, 4)) +
  labs(title = 'estimated GP')
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
ggsave('ypredicted_ordinal.png', width = 6, height = 4)

