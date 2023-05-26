functions {
  //Square root vector of spectral densities (for a squared exponential kernel)
  vector sqrt_diagSPD_nD(real gpscale, real lscale, vector L, matrix indices, int D) {
    return gpscale *  sqrt(sqrt(2*pi())^D * (lscale)) * exp(-.25 * (indices^2 * (lscale*pi() ./ (2*L))^2));
  }
	
	//Matrix of eigenfunction values
  matrix PHI_2D(int N, int M1, int M2, real L1, real L2, vector x1, vector x2) {
    matrix[N,M1*M2] PHI;
    matrix[N,M1] PHI_1 = sin(diag_post_multiply(rep_matrix(pi()/(2*L1) * (x1+L1), M1), linspaced_vector(M1, 1, M1)))/sqrt(L1);
    matrix[N,M2] PHI_2 = sin(diag_post_multiply(rep_matrix(pi()/(2*L2) * (x2+L2), M2), linspaced_vector(M2, 1, M2)))/sqrt(L2);
    PHI[,1:M2] = rep_matrix(PHI_1[,1], M2) .* PHI_2;
    for(i in 2:M1)
      PHI[,1:(M2*i)] = append_col(PHI[,1:(M2*(i-1))], rep_matrix(PHI_1[,i], M2) .* PHI_2);
    return PHI;
  }
  
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
}

data {
	int<lower=1> D; // dimensions
	vector[D] L; // boundary domain
	array[D] int<lower=1> M; // no. basis functions
	int<lower=1> M_nD; // n_M1 * n_M2 
	matrix[M_nD,D] indices; // basis functions for each dimension
	int<lower=1> N_sample;
	int N_pred;
	int vv_sample[N_sample];
	matrix[N_pred,D] x_pred; // coordinates
	int<lower=1> K; // no. ordinal categories
	array[N_pred] int<lower=1, upper=K> y_pred; // response
	
}

transformed data {
	matrix[N_pred,M_nD] PHI;
	PHI = PHI_2D(N_pred, M[1], M[2], L[1], L[2], x_pred[,1], x_pred[,2]);
}

parameters {
	vector[M_nD] beta; // zscore for GP?
	real<lower=0> lscale; // GP length scale
	real<lower=0> gpscale; // GP magnitude
	ordered[K-1] cutpoint; // internal cut points
}

transformed parameters{
	vector[N_sample] f;
	vector[M_nD] SPD_beta;
 {
	vector[M_nD] SPD;
	SPD = sqrt_diagSPD_nD(gpscale, lscale, L, indices, D);
	SPD_beta = SPD .* beta;
	f= PHI[vv_sample,] * SPD_beta;
 }
}

model{
	beta ~ normal(0,1);
	lscale ~ inv_gamma(2,.5);//gamma(2, 3); // inv_gamma(2,.5);
	gpscale ~ normal(0,4);
	cutpoint ~ induced_dirichlet(rep_vector(1, K), 0);
	
	y_pred[vv_sample] ~ ordered_logistic(f, cutpoint);
}

generated quantities{
  vector[N_pred] f_pred;
  vector[N_sample] elpd;
  array[N_sample] int<lower=1, upper=K> y_rep;
  array[N_pred] int<lower=1, upper=K> y_rep_all;
  
  f_pred= PHI * SPD_beta;
	for(i in 1:N_sample){
	  elpd[i] = ordered_logistic_lpmf(y_pred[vv_sample][i] | f[i], cutpoint);
	} 
	  
	for(n in 1:N_sample) y_rep[n] = ordered_logistic_rng(f[n], cutpoint);
	for(n in 1:N_pred) y_rep_all[n] = ordered_logistic_rng(f_pred[n], cutpoint);
  
}


