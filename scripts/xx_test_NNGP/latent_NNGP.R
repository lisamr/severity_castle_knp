# can i get the latent model work in coordinates more similar to my data? 

library(tidyverse)
library(cmdstanr)

source('scripts/xx_test_NNGP/NNMatrix.R')


# define other functions --------------------------------------------------

rmvn <- function(N, mu = 0, V = matrix(1)){
  P <- length(mu)
  if(any(is.na(match(dim(V), P))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(N * P), ncol = P) %*% D + rep(mu, rep(N, P)))
}

make_data <- function(N, B, sigma.sq, tau.sq, phi, coordmax=1){
  coords <- cbind(runif(N, min = 0, max = coordmax), 
                  runif(N, min = 0, max = coordmax))
  X <- as.matrix(cbind(1, rnorm(N)))
  B <- as.matrix(c(1, B))
  D <- as.matrix(dist(coords))
  R <- exp(-phi*D)
  w <- rmvn(1, rep(0, N), sigma.sq*R)
  Y <- rnorm(N, X %*% B + w, sqrt(tau.sq))
  return(list(
    coords = coords, 
    X = X, 
    Y = Y
  ))
}

prep_data <- function(N, B, sigma.sq, tau.sq, phi, M, coordmax=1){
  
  # parameters defining priors
  P = 1                  # number of regression coefficients
  uB = rep(0, P + 1)     # mean vector in the Gaussian prior of beta
  VB = diag(P + 1)*1000  # covariance matrix in the Gaussian prior of beta
  ss = 3 * sqrt(sigma.sq)       # scale parameter in the normal prior of sigma 
  st = 3 * sqrt(tau.sq)     # scale parameter in the normal prior of tau     
  ap = 3; bp = 0.5       # shape and scale parameter in the inv-gamma prior of phi 
  
  # generate GP data 
  dat <- make_data(N, B, sigma.sq, tau.sq, phi, coordmax)
  
  # neighbor index
  NN.matrix <- NNMatrix(coords = dat$coords, n.neighbors = M, n.omp.threads = 2)
  
  # prep and run model 
  data_list <- list(N = N, M = M, P = P, 
                    # sorted Y and X
                    Y = dat$Y[NN.matrix$ord], 
                    X = dat$X[NN.matrix$ord, ], 
                    NN_ind = NN.matrix$NN_ind, 
                    NN_dist = NN.matrix$NN_dist, 
                    NN_distM = NN.matrix$NN_distM,  
                    uB = uB, VB = VB, sS = ss, sT = st, ap = ap, bp = bp)
  myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.5, phi = 12, 
                      w_b1 = rep(0, N)), 
                 list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5, 
                      w_b1 = rep(0.1, N)), 
                 list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9 ,
                      w_b1 = rep(0, N)),
                 list(beta = c(2, 2), sigma = .5, tau = 0.3, phi = 4 ,
                      w_b1 = rep(0, N)))
  
  return(list(data_list = data_list, myinits = myinits))
}





fit_models <- function(N, B, sigma.sq, tau.sq, phi, M, n_chains, iters, coordmax=1){
  
  prepped_data <- prep_data(N, B, sigma.sq, tau.sq, phi, M, coordmax)
  
  t1 <- proc.time()
  fit_b1 <- nngp_latentb1$sample(prepped_data$data_list, init = prepped_data$myinits,
                                 iter_sampling = iters, iter_warmup = iters, 
                                 parallel_chains = n_chains, chains = n_chains,
                                 output_dir = 'scripts/xx_test_NNGP/model_fits/',
                                 output_basename = glue::glue("N{N}_M{M}")
  )
  t2 <- proc.time()
  tdelta <- (t2 - t1)[3]
  
  return(list(fit_b1, tdelta = tdelta))
}



# set parameters ----------------------------------------------------------

if(!dir.exists('scripts/xx_test_NNGP/model_fits')) dir.create('scripts/xx_test_NNGP/model_fits')


set.seed(1234)
# data generating parameters
N <- 500 #c(100, 200, 500, 1000)
B <- 5
sigma.sq <- 2
tau.sq <- 0.1
phi <- 3
M = 6
coordmax <- 10

# look at the kernel
exponential <- function(d, etasq, rho) etasq * exp(-(d * rho))
curve(exponential(x, sigma.sq, 3), 0, 10)
curve(exponential(x, sigma.sq, 3), 0, 1)

prepped_data <- prep_data(N, B, sigma.sq, tau.sq, phi, M, 10)
sort(prepped_data$data_list$NN_dist[50,])

# parameters for sampling
n_chains = 4
iters = 500

# load up model
nngp_latentb1 <- cmdstan_model('scripts/xx_test_NNGP/latent_normal_tutorial.stan')

results <- fit_models(N, B, sigma.sq, tau.sq, phi, M, n_chains, iters,
           coordmax)



# look at model summaries -------------------------------------------------

# do they estimate parameters ok at least? kind of. tau is off. might need ot increase M. 

csv_files <- fs::dir_ls('scripts/xx_test_NNGP/model_fits', glob = '*N500*')
fit2 <- as_cmdstan_fit(csv_files)
fit2$summary(c('beta', 'sigmasq', 'tausq', 'phi'))

# variable  mean median     sd    mad      q5   q95  rhat ess_bulk
# <chr>    <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>    <dbl>
# 1 beta[1]  0.298  0.340 0.749  0.612  -0.910  1.37   1.01    1685.
# 2 beta[2]  5.00   5.00  0.0296 0.0296  4.95   5.05   1.00     797.
# 3 sigmasq  3.04   2.62  1.67   0.896   1.65   5.90   1.00     819.
# 4 tausq    0.133  0.131 0.0286 0.0282  0.0880 0.184  1.02     123.
# 5 phi      3.60   3.54  1.37   1.45    1.48   5.89   1.00     579.
