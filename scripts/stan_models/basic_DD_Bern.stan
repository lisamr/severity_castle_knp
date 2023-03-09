data{
  int N; // number of sites
  int M; // number of rings. M > N
  matrix[N,M] dd_var; // mean of dist-dep var in each ring
  vector<lower=0>[M] area; // area of rings. Rings same for all sites, so vector of area and distance applies to all sites.
  vector<lower=0>[M] dist_sq; // distance^2 from ring m to site centers
  array[N] int<lower=0,upper=1> y; // response variable
  int Nsim;
  vector[Nsim] dist_sq_sim;
  vector[Nsim] area_sim;
}
parameters{
  real<lower=0> delta;
  real a0;
  real beta;
}
transformed parameters{
  vector<lower=0>[M] w0;
  vector<lower=0>[M] w;
  vector[N] DD;
  vector[N] p;
  
  // calculate weights parameter
  w0 = exp(-dist_sq / (2 * square(delta))) .* area; // .* allows for row-wise mulitplication
  w = w0 / sum(w0);
  
  // distance dependent effect = sum of each ring's effect
  DD = dd_var * w; // matrix multiplication ~ for(i in 1:N) sum(dd_var[i,] * w)
  
  // main model
  p = a0 + beta * DD;
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  delta ~ lognormal(0, 1);
  
  y ~ bernoulli_logit(p);
}
generated quantities{
  array[N] int y_rep;
  vector[Nsim] w_rep;
  for(i in 1:N){
    y_rep[i] = bernoulli_logit_rng(p[i]);
  }
  
  // generate weights for simulated data
  w_rep = exp(-dist_sq_sim / (2 * square(delta))) .* area_sim; 
  w_rep = w_rep / sum(w_rep);
}
