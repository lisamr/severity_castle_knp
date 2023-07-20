// ordinal regression + distance-dependent predictors

functions {
  // from Betancourt's page on ordinal regressions
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
  // dd-var stuff
  int nDD; // number of DD variables
  int M; // number of rings. M > N
  array [nDD] matrix[N,M] dd_var; // mean of dist-dep var in each ring
  vector<lower=0>[M] area; // area of rings. Rings same for all sites, so vector of area and distance applies to all sites.
  vector<lower=0>[M] dist_sq; // distance^2 from ring m to site centers
  int<lower=1> nInt;
  array[nInt] int<lower=1> int_var;
  // // for simulating posterior predictions
  int<lower=1> Nsim;
  vector[Nsim] Xsim;
  matrix[Nsim,M] dd_var_sim;
}
transformed data{
  real L = 0.01; // lower limit for delta (10m). Prevents w from being NaN
}
parameters {
  vector[J] B;
  vector[nDD] beta;
  array[nInt] vector[nDD] beta_int;
  ordered[K - 1] c; // cut points, i.e. log-cumulative-probability intercepts
  vector<lower=L>[nDD] delta;
}
transformed parameters{
  vector[N] phi;
  array[nDD] vector<lower=0>[M] w0;
  array[nDD] vector<lower=0>[M] w;
  matrix[N,nDD] DD;
  vector[N] Beta_int;
  
  for(i in 1:nDD){
    // calculate weights parameter
    w0[i] = exp(-dist_sq / (2 * square(delta[i]))) .* area; 
    w[i] = w0[i] / sum(w0[i]);
    // distance dependent effect = sum of each ring's effect
    DD[,i] = dd_var[i] * w[i];
  }
  
  // Calculate sum of interaction effects
  Beta_int = rep_vector(0, N);
  for (i in 1:nInt) {
    matrix[N,nDD] interaction_term = DD .* rep_matrix(X[,int_var[i]], nDD);
    Beta_int = Beta_int + interaction_term * beta_int[i];
  }
  

  // main model
  phi = X*B + DD*beta + Beta_int;
}
model {
  
  // Prior model
  B ~ normal(0, 1);
  beta ~ normal(0, 1);
  c ~ induced_dirichlet(rep_vector(1, K), 0); // parameter alpha, second argument is anchor point
  for(i in 1:nDD) delta[i] ~ lognormal(-1, 1.5)T[L,];
  
  // Observational model
  y ~ ordered_logistic(phi, c);
}
generated quantities {
  array[N] int<lower=1, upper=K> y_rep; 
  vector[Nsim] phi_sim;
  vector[Nsim] DD_sim;
  for(n in 1:N) y_rep[n] = ordered_logistic_rng(phi[n], c);
  
  // estimate log-cumulative-odds of getting each category for the sim data.
  // simulated distance dependent effect (assumes first effect)
  DD_sim = dd_var_sim * w[1];
  // assumes all variables constant except 1st interaction and 1st DD vars
  // beta_int[1][1] <-- effect of first int var on first DD var
  phi_sim = Xsim*B[int_var[1]] + DD_sim*beta[1] + (DD_sim.*Xsim)*beta_int[1][1];
  

}