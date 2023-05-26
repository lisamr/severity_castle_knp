functions{
   // for ordinal response variable
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
  
  vector gp_pred_rng(vector[] x_pred,
                     vector y1, vector[] x,
                     real magnitude, real[] length_scale) {
    int N = rows(y1);
    int N_pred = size(x_pred);
    vector[N_pred] f2;
    {
      matrix[N, N] K = gp_matern52_cov(x, magnitude, length_scale);
      matrix[N, N] L_K = cholesky_decompose(K);
      vector[N] L_K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N] K_div_y1 = mdivide_right_tri_low(L_K_div_y1', L_K)';
      matrix[N, N_pred] k_x_x_pred = gp_matern52_cov(x, x_pred, magnitude, length_scale);
      f2 = (k_x_x_pred' * K_div_y1);
    }
    return f2;
  }
  
}
data{
	int<lower=1> N_sample;
	int N_pred;
	array[N_sample] int idx;
  int<lower=1> nvars;
  int<lower=1> K; // no. ordinal categories
	array[N_pred] int<lower=1, upper=K> y; // response
  array[N_pred] vector[2] coords;
  matrix[N_pred,nvars] X;
  // dd-var stuff
  int nDD; // number of DD variables
  int R; // number of rings. M > N
  array [nDD] matrix[N_pred,R] dd_var; // mean of dist-dep var in each ring
  vector<lower=0>[R] area; // area of rings. Rings same for all sites, so vector of area and distance applies to all sites.
  vector<lower=0>[R] dist_sq; // distance^2 from ring m to site centers

}
transformed data{
  real L = 0.01; // lower limit for delta (10m). Prevents w from being NaN
}
parameters{
  vector[nvars] B;
  ordered[K-1] cutpoint; 
  //GP terms
  vector[N_sample] zGP;
  real<lower=0> eta;
  real<lower=0> rho;
  // dist-dep stuff
  vector[nDD] beta;
  ordered[K - 1] c; // cut points, i.e. log-cumulative-probability intercepts
  vector<lower=L>[nDD] delta;
}
transformed parameters{ 
  vector[N_sample] phi; 
  vector[N_sample] GP;
  array[nDD] vector<lower=0>[R] w0;
  array[nDD] vector<lower=0>[R] w;
  matrix[N_sample,nDD] DD;
  
  for(i in 1:nDD){
    // calculate weights parameter
    w0[i] = exp(-dist_sq / (2 * square(delta[i]))) .* area; 
    w[i] = w0[i] / sum(w0[i]);
    // distance dependent effect = sum of each ring's effect
    DD[,i] = dd_var[i][idx,] * w[i];
  }
  
  // Define the GP matrices
  {
    matrix[N_sample, N_sample] gp_K;
    matrix[N_sample, N_sample] gp_KL;
    
    gp_K = gp_exp_quad_cov(coords[idx], eta, rho);
    gp_KL = cholesky_decompose(add_diag(gp_K, 1e-9));
    GP = gp_KL * zGP;
  } 
  
  //main model
  phi = X[idx,]*B + DD*beta + GP; 

}
model{
  B ~ normal(0, 2);
  cutpoint ~ induced_dirichlet(rep_vector(1, K), 0);
  beta ~ normal(0, 2);
  for(i in 1:nDD) delta[i] ~ lognormal(-1, 1.5)T[L,]; //could replace with inv-gamma
  zGP ~ std_normal();
  rho ~ inv_gamma(2, .5); 
  eta ~ std_normal();
	
	y[idx] ~ ordered_logistic(phi, cutpoint);
}
generated quantities{
  vector[N_sample] elpd;
  array[N_sample] int<lower=1, upper=K> y_rep;
  
	for(i in 1:N_sample) elpd[i] = ordered_logistic_lpmf(y[idx][i] | phi[i], cutpoint);
	for(i in 1:N_sample) y_rep[i] = ordered_logistic_rng(phi[i], cutpoint);
}
