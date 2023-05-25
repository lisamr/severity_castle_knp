# build ordinal regression with NNGP

rm(list = ls())
library(tidyverse)
library(cmdstanr)

source('scripts/xx_test_NNGP/NNMatrix.R')




# functions ---------------------------------------------------------------

# kernel function
exp_quad <- function(d, eta, rho) eta^2 * exp(-(d^2) / (2 * rho^2))
exponential <- function(d, eta, rho) eta^2 * exp(-(d * rho))

make_data <- function(N, eta, sd, rho, B, kernel, coordmax = 10){
  coords <- cbind(runif(N, 0, coordmax), runif(N, 0, coordmax))
  X <- cbind(1, rnorm(N))
  
  # deal with the GP
  dmat <- fields::rdist(coords)
  K <- do.call(kernel, list(dmat, eta, rho))
  # cholesky decomposition
  zGP <- matrix(rnorm(N), nrow = N, ncol = 1)
  diag(K) <- diag(K) + 1e-9
  GP <- t(chol(K)) %*% zGP
  
  Y <- rnorm(N, X %*% B + GP, sd)
  return(list(
    coords = coords, 
    X = X, 
    Y = Y
  ))
}



# simulate data -----------------------------------------------------------

# first work with normal data
N <- 500 
eta <- 1
sigma <- .5
rho <- 3/.5
beta <- c(1, 0)
M <- 10
curve(exponential(x, eta, rho), 0, 1)
simdata <- make_data(N = N, eta, sigma, rho, B = beta, kernel = 'exponential', coordmax = 1)
dat_tibble <- with(simdata, 
            tibble(
              x = coords[,1],
              y = coords[,2],
              Y = Y)
            )


ggplot(dat_tibble, aes(x, y)) +
  geom_point(aes(color = Y, size = Y), alpha = .9) +
  coord_equal()


# neighbor index
NN.matrix <- NNMatrix(coords = simdata$coords, n.neighbors = M, 
                      n.omp.threads = 2)
Check_Neighbors(NN.matrix$coords.ord, M, NN.matrix, 200)




# # tune inverse gamma prior
# stan_tune <- cmdstan_model('scripts/xx_test_NNGP/invG_tune.stan')
# drange <- range(NN.matrix$NN_dist[NN.matrix$NN_dist != 0])
# stan_tune$sample(list(l =drange[1], u = drange[2]), iter_warmup = 0, iter_sampling = 1, chains = 1, fixed_param = T)
# stan_tune$sample(list(l =2, u = 20), iter_warmup = 0, iter_sampling = 1, chains = 1, fixed_param = T)



# fit data ----------------------------------------------------------------

# model from the tutorial
stan_response <- cmdstan_model('scripts/xx_test_NNGP/response_NNGP.stan')

# shape and scale parameter in the inv-gamma prior of phi 
P = 1                  # number of regression coefficients
uB = rep(0, P + 1)     # mean vector in the Gaussian prior of beta
VB = diag(P + 1)*1000  # covariance matrix in the Gaussian prior of beta
ss = 3 * sqrt(eta)       # scale parameter in the normal prior of sigma 
st = 3 * sqrt(sigma)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       

data_list <- list(N = N, 
                  M = M, 
                  P = P,
                  #J = length(beta), 
                  # sorted Y and X
                  Y = simdata$Y[NN.matrix$ord],
                  X = cbind(simdata$X[NN.matrix$ord,]), 
                  NN_ind = NN.matrix$NN_ind, 
                  NN_dist = NN.matrix$NN_dist, 
                  NN_distM = NN.matrix$NN_distM,  
                  uB = uB, VB = VB, ss = ss, st = st,
                  ap = ap, bp = bp)

myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.5, phi = 12), 
               list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5), 
               list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9),
               list(beta = c(1, 1), sigma = 1, tau = 0.1, phi = 3))
# 
# myinits <-list(list(beta = c(1, 5), eta = 1, sigma = 0.5, rho = 12, 
#                     w = rep(.3, N), a0 = 1), 
#                list(beta = c(5, 5), eta = 1.5, sigma = 0.2, rho = 5, 
#                     w = rep(0.1, N), a0 = 2), 
#                list(beta = c(0, 0), eta = 2.5, sigma = 0.1, rho = 9 ,
#                     w = rep(0, N), a0 = -1),
#                list(beta = c(2, 2), eta = 1, sigma = 0.2, rho = 2, 
#                     w = rep(0.1, N), a0 = -2))

# N <- 500 
# eta <- 1
# sigma <- .5
# rho <- 6
# beta <- 0

fit_response <- stan_response$sample(data_list, init = myinits, 
                                 parallel_chains = 4, 
                                 iter_sampling = 500, iter_warmup = 500)
fit_response$summary(c('beta', 'sigma', 'tau', 'phi'))



# alter model: beta isn't multivariate ------------------------------------

# beta isn't multivariate. and change the names of the parametesr 
stan_normal <- cmdstan_model('scripts/xx_test_NNGP/marginal_normal.stan')

# shape and scale parameter in the inv-gamma prior of phi 
P = 1                  # number of regression coefficients
se = 3 * sqrt(eta)       # scale parameter in the normal prior of sigma 
ss = 3 * sqrt(sigma)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       

data_list <- list(N = N, 
                  M = M, 
                  P = P,
                  #J = length(beta), 
                  # sorted Y and X
                  Y = simdata$Y[NN.matrix$ord],
                  X = cbind(simdata$X[NN.matrix$ord,]), 
                  NN_ind = NN.matrix$NN_ind, 
                  NN_dist = NN.matrix$NN_dist, 
                  NN_distM = NN.matrix$NN_distM,  
                  se = se, ss = ss,
                  ap = ap, bp = bp,
                  # for matching with last data list
                  uB = uB, VB = VB, ss = ss, st = st)


fit_normal <- stan_normal$sample(data_list, 
                                     parallel_chains = 4, 
                                     iter_sampling = 500, iter_warmup = 500)
fit_normal$summary(c('beta', 'sigma', 'eta', 'rho'))
bayesplot::mcmc_pairs(fit_normal$draws(c('beta', 'sigma', 'eta', 'rho')))



# change the coordinates of the data --------------------------------------

# assume distance in units of km
N <- 500 
eta <- 1
sigma <- .5
rho <- 3
beta <- c(1, 0)
M <- 10
# look at kernel and prior for rho
curve(exponential(x, eta, rho), 0, 10)

simdata <- make_data(N = N, eta, sigma, rho, B = beta, 
                     kernel = 'exponential', coordmax = 10)
NN.matrix <- NNMatrix(coords = simdata$coords, n.neighbors = M, 
                      n.omp.threads = 2)

# shape and scale parameter in the inv-gamma prior of phi 
P = 1                  # number of regression coefficients
se = 3 * sqrt(eta)       # scale parameter in the normal prior of sigma 
ss = 3 * sqrt(sigma)     # scale parameter in the normal prior of tau     
ap = 2; bp = 1.5       

data_list <- list(N = N, 
                  M = M, 
                  P = P,
                  #J = length(beta), 
                  # sorted Y and X
                  Y = simdata$Y[NN.matrix$ord],
                  X = cbind(simdata$X[NN.matrix$ord,]), 
                  NN_ind = NN.matrix$NN_ind, 
                  NN_dist = NN.matrix$NN_dist, 
                  NN_distM = NN.matrix$NN_distM,  
                  se = se, ss = ss,
                  ap = ap, bp = bp,
                  # for matching with last data list and tutoral code
                  uB = uB, VB = VB, sS = se, sT = ss)

fit_normal_bigger <- stan_normal$sample(data_list, 
                                 parallel_chains = 4, 
                                 iter_sampling = 500, iter_warmup = 500)
fit_normal_bigger$summary(c('beta', 'sigma', 'eta', 'rho'))
#bayesplot::mcmc_pairs(fit_normal_bigger$draws(c('beta', 'sigma', 'eta', 'rho')))




# try out inv-gamma prior on decay par ------------------------------------

# tune inverse gamma prior
stan_tune <- cmdstan_model('scripts/xx_test_NNGP/invG_tune.stan')
stan_tune$sample(list(l =1, u = 6), iter_warmup = 0, iter_sampling = 1, chains = 1, fixed_param = T)
curve(pscl::densigamma(x, 7.3, 15), 0.01, 10) # solid = inv-gamma

data_list_invgamma <- data_list
data_list_invgamma$ap <- 7.3
data_list_invgamma$bp <- 15

stan_normal_invg <- cmdstan_model('scripts/xx_test_NNGP/marginal_normal_invgamma.stan')
fit_normal_invg <- stan_normal_invg$sample(data_list_invgamma, 
                                        parallel_chains = 4, 
                                        iter_sampling = 1000, iter_warmup = 1000)
fit_normal_invg$summary(c('beta', 'sigma', 'eta', 'rho'))




# do latent model ---------------------------------------------------------

# this isnt working...low ebfmi and bad sampling
stan_normal_latent <- cmdstan_model('scripts/xx_test_NNGP/latent_normal.stan')
fit_latent <- stan_normal_latent$sample(data_list_invgamma, 
                                        parallel_chains = 4, 
                                        iter_sampling = 500, iter_warmup = 500)
fit_latent$summary(c('beta', 'sigma', 'eta', 'rho'))
#4:17

# does it work with the one from the tutorial?
stan_normal_latent_tut <- cmdstan_model('scripts/xx_test_NNGP/latent_normal_tutorial.stan')
fit_latent_tut <- stan_normal_latent$sample(data_list, 
                                        parallel_chains = 4, 
                                        iter_sampling = 500, iter_warmup = 500)
fit_latent_tut$summary(c('beta', 'sigma', 'eta', 'rho'))


# ordinal model with exponential kernel -----------------------------------

# couldn't figure out how to use the gaussian kernel, so stick with exponential






# 
# # gaussian kernel ---------------------------------------------------------
# 
# # data needs to be resimulated
# # stan model needs to be modified at the kernel
# 
# N <- 500 
# eta <- 1
# sigma <- .5
# rho <- 2
# beta <- c(0, .5)
# M <- 10
# curve(exp_quad(x, eta, rho), 0, 10)
# simdata <- make_data(N = N, eta, sigma, rho, B = beta, kernel = 'exp_quad', 
#                      coordmax = 10)
# dat_tibble <- with(simdata, 
#                    tibble(
#                      x = coords[,1],
#                      y = coords[,2],
#                      Y = Y)
# )
# 
# 
# ggplot(dat_tibble, aes(x, y)) +
#   geom_point(aes(color = Y, size = Y), alpha = .9) +
#   coord_equal()
# 
# 
# # neighbor index
# NN.matrix <- NNMatrix(coords = simdata$coords, n.neighbors = M, 
#                       n.omp.threads = 2)
# Check_Neighbors(NN.matrix$coords.ord, M, NN.matrix, 200)
# 
# 
# # tune inverse gamma prior
# stan_tune <- cmdstan_model('scripts/xx_test_NNGP/invG_tune.stan')
# drange <- range(NN.matrix$NN_dist[NN.matrix$NN_dist != 0])
# stan_tune$sample(list(l =.7, u = drange[2]), iter_warmup = 0, iter_sampling = 1, chains = 1, fixed_param = T)
# curve(pscl::densigamma(x, 3.85, 6.87), 0.01, 10) # solid = inv-gamma
# curve(dgamma(x, 2, 1), .01, 10, add = T, lty = 2) # dashed = gamma
# 
# # pars for priors
# ss = 3 * sqrt(eta)       # scale parameter in the normal prior of sigma 
# st = 3 * sqrt(sigma)     # scale parameter in the normal prior of tau     
# ap = 2; bp = 1     
# 
# data_list <- list(N = N, 
#                   M = M, 
#                   P = length(beta)-1,
#                   #J = length(beta), 
#                   # sorted Y and X
#                   Y = simdata$Y[NN.matrix$ord],
#                   X = cbind(simdata$X[NN.matrix$ord,]), 
#                   NN_ind = NN.matrix$NN_ind, 
#                   NN_dist = NN.matrix$NN_dist, 
#                   NN_distM = NN.matrix$NN_distM,  
#                   ss = ss, st = st,
#                   ap = ap, bp = bp)
# 
# # fit again with exp_quad kernel
# stan_normal_expquad <- cmdstan_model('scripts/xx_test_NNGP/marginal_normal_expquad.stan')
# 
# fit_normal_expquad <- stan_normal_expquad$sample(data_list, init = myinits, 
#                                                  parallel_chains = 4, 
#                                                  iter_sampling = 500, iter_warmup = 500)
# 
# # not looking good. intercept can't be estimated, not good overall
# fit_normal_expquad$summary(c('beta', 'sigma', 'tau', 'phi'))
# bayesplot::mcmc_pairs(fit_normal_expquad$draws(c('beta', 'sigma', 'tau', 'phi')))
# 
# 
