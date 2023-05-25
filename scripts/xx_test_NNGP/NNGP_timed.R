# follow along Lu Zhang's script
# need to use latent model centered at intercept to accomodate ordinal response variable

# need to assess time burden...loop over multiple sample sizes and compare time.



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

make_data <- function(N, B, sigma.sq, tau.sq, phi){
  coords <- cbind(runif(N), runif(N))
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


fit_models <- function(N, B, sigma.sq, tau.sq, phi, M, n_chains, iters){
  
  # parameters defining priors
  P = 1                  # number of regression coefficients
  uB = rep(0, P + 1)     # mean vector in the Gaussian prior of beta
  VB = diag(P + 1)*1000  # covariance matrix in the Gaussian prior of beta
  ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
  st = 3 * sqrt(0.1)     # scale parameter in the normal prior of tau     
  ap = 3; bp = 0.5       # shape and scale parameter in the inv-gamma prior of phi 
  
  # generate GP data 
  dat <- make_data(N, B, sigma.sq, tau.sq, phi)
  
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
  
  
  t1 <- proc.time()
  fit_b1 <- nngp_latentb1$sample(data_list, init = myinits,
                                 iter_sampling = iters, iter_warmup = iters, 
                                 parallel_chains = n_chains, chains = n_chains,
                                 output_dir = 'scripts/xx_test_NNGP/model_fits/',
                                 output_basename = glue::glue("N{N}_M{M}")
  )
  t2 <- proc.time()
  tdelta <- (t2 - t1)[3]
  
  return(list(N = N, M = M, tdelta = tdelta))
}
 


# set parameters ----------------------------------------------------------

if(!dir.exists('scripts/xx_test_NNGP/model_fits')) dir.create('scripts/xx_test_NNGP/model_fits')

set.seed(1234)
# data generating parameters
N <- 500 #c(100, 200, 500, 1000)
B <- 5
sigma.sq <- 2
tau.sq <- 0.1
phi <- 3 / 0.5
M = 6

# parameters for sampling
n_chains = 4
iters = 500

# load up model
nngp_latentb1 <- cmdstan_model('scripts/xx_test_NNGP/latent_normal_tutorial.stan')


results <- map(N, ~ fit_models(.x, B, sigma.sq, tau.sq, phi, M, n_chains, iters))

results_df <- transpose(results) %>% 
  map(unlist) %>% 
  map(cbind) %>% 
  do.call(cbind, .) %>% 
  as_tibble() %>% 
  set_names(c('N', 'M', 'tdelta'))

ggplot(results_df, aes(N, log(tdelta))) +
  geom_point() +
  geom_line()

mod <- lm(log(tdelta) ~ N, data = results_df)
v <- c(N, 1000, 2000, 5000)
secs <- predict(mod, newdata = data.frame(N = v)) %>% exp
hrs <- secs/60/60
plot(v, log(hrs), type = 'o')

# look at model summaries -------------------------------------------------

# do they estimate parameters ok at least? kind of. tau is off. might need ot increase M. 
fit <- cmdstanr::as_cmdstan_fit(c('scripts/xx_test_NNGP/model_fits/N1000_M6-1.csv',
                           'scripts/xx_test_NNGP/model_fits/N1000_M6-2.csv',
                           'scripts/xx_test_NNGP/model_fits/N1000_M6-3.csv',
                           'scripts/xx_test_NNGP/model_fits/N1000_M6-4.csv'))
fit$summary(c('beta', 'sigma', 'tau', 'phi'))


csv_files <- fs::dir_ls('scripts/xx_test_NNGP/model_fits', glob = '*N500*')
fit2 <- as_cmdstan_fit(csv_files)
fit2$summary(c('beta', 'sigmasq', 'tausq', 'phi'))


