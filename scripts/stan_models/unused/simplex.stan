data{
  int N; // samples
  int M; // rings. M > N
  matrix[N,M] dd_var; // mean of d-d var in each ring
  vector[N] y; // response variable
}
transformed data{
  vector[M] w_prior = rep_vector(2.0, M);
}
parameters{
  simplex[M] w; // must sum to 1? contrain to be positive.
  real a0;
  real beta;
  real<lower=0> sigma;
}
transformed parameters{
  vector[N] DD;
  vector[N] mu;

  // distance dependent effect = sum of each ring's effect
  DD = dd_var * w;
  
  // main model
  mu = a0 + beta * DD;
  
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  w ~ dirichlet(w_prior);
  
  y ~ normal(mu, sigma);
}
generated quantities{
  vector[N] y_rep;
  for(i in 1:N){
    y_rep[i] = normal_rng(mu[i], sigma);
  }
}
