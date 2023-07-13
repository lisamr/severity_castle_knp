functions {
  //Square root vector of spectral densities (for a squared exponential kernel)
  vector sqrt_diagSPD_nD(real gpscale, vector lscale, vector L, matrix indices, int D) {
    return gpscale *  sqrt(sqrt(2*pi())^D * prod(lscale)) * exp(-.25 * (indices^2 * (lscale*pi() ./ (2*L))^2));
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
	int<lower=1> N;
	int N_total;
	int idx[N];
	matrix[N_total,D] coords; // coordinates
	int<lower=1> J; // number of predictors
  matrix[N_total,J] X; // data matrix of predictors
	int<lower=1> K; // no. ordinal categories
	array[N_total] int<lower=1, upper=K> y; // response
	
}

transformed data {
	matrix[N_total,M_nD] PHI;
	PHI = PHI_2D(N_total, M[1], M[2], L[1], L[2], coords[,1], coords[,2]);
}

parameters {
	vector[M_nD] zGP; // zscore for GP?
	vector<lower=0>[D] lscale; // GP length scale
	real<lower=0> gpscale; // GP magnitude
	vector[J] B;
	ordered[K-1] cutpoint; // internal cut points
}

transformed parameters{
	vector[N] f;
	vector[M_nD] SPD_GP;
	vector[N] mu;
 {
	vector[M_nD] SPD;
	SPD = sqrt_diagSPD_nD(gpscale, lscale, L, indices, D);
	SPD_GP = SPD .* zGP;
	f= PHI[idx,] * SPD_GP;
 }
 
 mu = f + X[idx,]*B;
}

model{
	zGP ~ normal(0,1);
	lscale ~ inv_gamma(2,.5);//gamma(2, 3); // inv_gamma(2,.5);
	gpscale ~ normal(0,4);
	B ~ normal(0, 2);
	cutpoint ~ induced_dirichlet(rep_vector(1, K), 0);
	
	y[idx] ~ ordered_logistic(mu, cutpoint);
}

generated quantities{
  vector[N_total] mu_pred;
  vector[N_total] f_pred;
  vector[N] log_lik;
  array[N] int<lower=1, upper=K> y_rep;
  array[N_total] int<lower=1, upper=K> y_rep_all;
  
  f_pred = PHI * SPD_GP;
  mu_pred = f_pred + X*B;
  
	for(i in 1:N){
	  log_lik[i] = ordered_logistic_lpmf(y[idx][i] | mu[i], cutpoint);
	} 
	  
	for(n in 1:N) y_rep[n] = ordered_logistic_rng(mu[n], cutpoint);
	for(n in 1:N_total) y_rep_all[n] = ordered_logistic_rng(mu_pred[n], cutpoint);
  
}


