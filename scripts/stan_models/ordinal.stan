// from Betancourt's page on ordinal regressions
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
data {
  int<lower=1> N; // Number of observations
  int<lower=1> K; // Number of ordinal categories
  int<lower=1> J; // number of predictors
  matrix[N,J] X; // data matrix of predictors
  array[N] int<lower=1, upper=K> y; // Observed ordinals
  // for simulating posterior predictions
  int<lower=1> Nsim; 
  matrix[Nsim,J] Xsim; 
}
/*transformed data{
  real gamma = 0;
}*/
parameters {
  vector[J] B;
  ordered[K - 1] c; // (Internal) cut points
}
transformed parameters{
  vector[N] phi;
  phi = X*B;
}
model {
  // Prior model
  B ~ normal(0, 1);
  c ~ induced_dirichlet(rep_vector(1, K), 0); // parameter alpha, second argument is anchor point

  // Observational model
  y ~ ordered_logistic(phi, c);
}
generated quantities {
  array[N] int<lower=1, upper=K> y_rep;
  vector[Nsim] phi_sim;
  for(n in 1:N) y_rep[n] = ordered_logistic_rng(phi[n], c);
  
  // estimate log-cumulative-odds of getting each category for the sim data.
  phi_sim = Xsim*B;
  
}
