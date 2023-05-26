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
	array[N_sample] int vv_sample;
  //int<lower=1> nvars;
  int<lower=1> K; // no. ordinal categories
	array[N_pred] int<lower=1, upper=K> y_pred; // response
  array[N_pred] vector[2] coords;
  //matrix[N_pred,nvars] X;
}
parameters{
  //vector[nvars] B;
  //GP terms
  vector[N_sample] zGP;
  real<lower=0> eta;
  real<lower=0> rho;
  ordered[K-1] cutpoint; 
}
transformed parameters{ 
  vector[N_sample] mu; 
  vector[N_sample] GP;
  
  // Define the GP matrices
  {
    matrix[N_sample, N_sample] gp_K;
    matrix[N_sample, N_sample] gp_KL;
    
    gp_K = gp_exp_quad_cov(coords[vv_sample], eta, rho);
    gp_KL = cholesky_decompose(add_diag(gp_K, 1e-9));
    GP = gp_KL * zGP;
  } 
  
  //main model
  mu = GP; //X[vv_sample,]*B + GP; 

}
model{
  //B ~ normal(0, 2);
  zGP ~ std_normal();
  rho ~ inv_gamma(2, .5); 
  eta ~ std_normal();
  cutpoint ~ induced_dirichlet(rep_vector(1, K), 0);
	
	y_pred[vv_sample] ~ ordered_logistic(mu, cutpoint);
}
generated quantities{
  vector[N_sample] elpd;
  array[N_sample] int<lower=1, upper=K> y_rep;
  //array[N_pred] int<lower=1, upper=K> y_rep_all;
  //vector[N_pred] gp_pred = gp_pred_rng(coords, y_pred[vv_sample], coords[vv_sample], eta, rho);
	for(i in 1:N_sample) elpd[i] = ordered_logistic_lpmf(y_pred[vv_sample][i] | mu[i], cutpoint);
	for(i in 1:N_sample) y_rep[i] = ordered_logistic_rng(mu[i], cutpoint);
	//for(n in 1:N_pred) y_rep_all[n] = ordered_logistic_rng(gp_pred[n], cutpoint);
}
