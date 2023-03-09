data{
  int N; // number of sites
  int M; // number of rings. M > N
  matrix[N,M] dd_var; // mean of dist-dep var in each ring
  vector<lower=0>[M-1] area; // area of rings. Rings same for all sites, so vector of area and distance applies to all sites.
  vector<lower=0>[M-1] dist_sq; // distance^2 from ring m to site centers
  vector[N] y; // response variable
}
parameters{
  real<lower=0> delta;
  real a0;
  real beta;
  real<lower=0> sigma;
}
transformed parameters{
  vector<lower=0>[M-1] w0;
  vector<lower=0>[M] w;
  vector[N] DD;
  vector[N] mu;
  
  // calculate weights parameter
  w0 = exp(-dist_sq / (2 * square(delta))) .* area; // .* allows for row-wise mulitplication
  w0 = w0 / sum(w0);
  w = append_row(0.0, w0);
  
  // distance dependent effect = sum of each ring's effect
  DD = dd_var * w; // matrix multiplication ~ for(i in 1:N) sum(dd_var[i,] * w)
  
  // main model
  mu = a0 + beta * DD;
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  delta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  y ~ normal(mu, sigma);
}
generated quantities{
  vector[N] y_rep;
  for(i in 1:N){
    y_rep[i] = normal_rng(mu[i], sigma);
  }
}
