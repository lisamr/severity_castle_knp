# follow along Lu Zhang's script
# need to use latent model to accomodate ordinal response variable
# Let N(si) be at most M closest points to S_i among the locations indexed 
# less than i 

library(tidyverse)
library(cmdstanr)

source('scripts/xx_test_NNGP/NNMatrix.R')


# quick visualization of M, N, index --------------------------------------


N = 1000
M = 100
pts <- tibble(x = runif(N, 0, 10), y = runif(N, 0, 10)) %>% 
  mutate(j = rank(x)) %>% 
  arrange(j)

ggplot(pts, aes(x, y, color = j)) +
  geom_point()

# find M nearest neighbor for a random point
dmat <- dist(cbind(pts$x, pts$y), upper = T) %>% as.matrix()



visNN <- function(i){
  
  neighbors <- which(pts$j < pts$j[i])
  mNN <- names(sort(dmat[i,neighbors])[1:M]) %>% as.integer
  
  category <- rep(0, N)
  category[mNN] <- 2
  category[i] <- 1
  pts$category <- as.factor(category)
  
  ggplot(pts, aes(x, y, color = category)) +
    geom_point() +
    coord_equal()
}

map(rev(500:510), visNN) 



# generate GP data --------------------------------------------------------

rmvn <- function(N, mu = 0, V = matrix(1)){
  P <- length(mu)
  if(any(is.na(match(dim(V), P))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(N * P), ncol = P) %*% D + rep(mu, rep(N, P)))
}

set.seed(1234)
N <- 500
coords <- cbind(runif(N), runif(N))
X <- as.matrix(cbind(1, rnorm(N)))
B <- as.matrix(c(1, 5))
sigma.sq <- 1
tau.sq <- .1
phi <- 3 / 0.5
#X = cbind(rep(1, N)) ; B = 1

D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0, N), sigma.sq*R)
Y <- rnorm(N, X %*% B + w, sqrt(tau.sq))

dat <- tibble(x = coords[,1],
       y = coords[,2],
       Y = Y, 
       w = w[,1])

ggplot(dat, aes(x, y)) +
 geom_point(aes(color = Y, size = w), alpha = .9) +
 coord_equal()
ggplot(dat, aes(w, Y)) +
  geom_point()



# build neighbor index ----------------------------------------------------


M = 6
NN.matrix <- NNMatrix(coords = coords, n.neighbors = M, n.omp.threads = 2, cov.model = 'exponential')
str(NN.matrix)
Check_Neighbors(NN.matrix$coords.ord, M, NN.matrix, 301)




# prep and run model ------------------------------------------------------

# parameters defining priors
P = 1                  # number of regression coefficients
uB = rep(0, P + 1)     # mean vector in the Gaussian prior of beta
VB = diag(P + 1)*1000  # covariance matrix in the Gaussian prior of beta
ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
st = 3 * sqrt(0.1)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       # shape and scale parameter in the inv-gamma prior of phi 




data_list <- list(N = N, M = M, P = P, 
             Y = Y[NN.matrix$ord], X = X[NN.matrix$ord, ],    # sorted Y and X
             NN_ind = NN.matrix$NN_ind, 
             NN_dist = NN.matrix$NN_dist, 
             NN_distM = NN.matrix$NN_distM,  
             uB = uB, VB = VB, ss = ss, st = st, ap = ap, bp = bp)

myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.5, phi = 12), 
               list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5), 
               list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9))

nngp_latent <- cmdstan_model('scripts/xx_test_NNGP/latent_poisson.stan')

fit <- nngp_latent$sample(data_list, 
                          iter_sampling = 300, iter_warmup = 300, 
                          parallel_chains = 4)



# try again recentering GP ------------------------------------------------


nngp_latentb1 <- cmdstan_model('scripts/xx_test_NNGP/latent_poisson_b1.stan')

myinits <-list(list(beta = c(1, 5), sigma = 1, tau = 0.5, phi = 12, 
                    w_b1 = rep(0, N)), 
               list(beta = c(5, 5), sigma = 1.5, tau = 0.2, phi = 5, 
                    w_b1 = rep(0.1, N)), 
               list(beta = c(0, 0), sigma = 2.5, tau = 0.1, phi = 9 ,
                    w_b1 = rep(0, N)))


fit_b1 <- nngp_latentb1$sample(data_list, init = myinits, 
                          iter_sampling = 300, iter_warmup = 300, 
                          parallel_chains = 3, chains = 3)

parameters <- c("beta", "sigmasq", "tausq", "phi", "w_b1")
fit_b1$summary(parameters)
