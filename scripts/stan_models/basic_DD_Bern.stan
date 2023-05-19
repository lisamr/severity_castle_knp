data{
  int N; // number of sites
  int M; // number of rings. M > N
  matrix[N,M] dd_var; // mean of dist-dep var in each ring
  vector<lower=0>[M] area; // area of rings. Rings same for all sites, so vector of area and distance applies to all sites.
  vector<lower=0>[M] dist_sq; // distance^2 from ring m to site centers
  array[N] int<lower=0,upper=1> y; // response variable
  int Nsim_c; // simulating counterfactual
  matrix[Nsim_c,M] dd_var_sim;
  int Nsim_d; // simulating distance dep functions
  vector[Nsim_d] dist_sq_sim;
  vector[Nsim_d] area_sim;
}
transformed data{
  real U = max(sqrt(dist_sq));
  real L = 0.03; //lower limit for delta (30m). prevents sampling issues. 
}
parameters{
  real<lower=L> delta; //real<lower=0, upper = U> delta;
  real a0;
  real beta;
  //real<lower=0> sigma;
}
transformed parameters{
  vector<lower=0>[M] w0;
  vector<lower=0>[M] w;
  vector[N] DD;
  vector[N] mu;
  
  // calculate weights parameter
  w0 = exp(-dist_sq / (2 * square(delta))) .* area; // .* allows for row-wise mulitplication
  w = w0 / sum(w0);
  
  if(is_nan(w[1])){
    print("uh-oh, w[1] isn't a number. w0[1] =", w0[1]);
  }
  
  // distance dependent effect = sum of each ring's effect
  DD = dd_var * w; // matrix multiplication ~ for(i in 1:N) sum(dd_var[i,] * w)
  
  // main model
  mu = a0 + beta * DD;
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  delta ~ lognormal(-1, 1.5)T[L,];
  //sigma ~ normal(0, 1);
  
  y ~ bernoulli_logit(mu);
}
generated quantities{
  array[N] int y_rep;
  for(i in 1:N){
    //y_rep[i] = normal_rng(mu[i], sigma);
    y_rep[i] = bernoulli_logit_rng(mu[i]);
  }
  
  // estimate counterfactual
  vector[Nsim_c] y_sim;
  for(i in 1:Nsim_c){
    vector[Nsim_c] DD_sim = dd_var_sim * w; 
    vector[Nsim_c] mu_sim = a0 + beta * DD_sim;
    y_sim[i] = bernoulli_logit_rng(mu_sim[i]);
  }
  
  
  // generate weights for simulated data
  vector[Nsim_d] w_rep;
  w_rep = exp(-dist_sq_sim / (2 * square(delta))) .* area_sim; 
  w_rep = w_rep / sum(w_rep);
}
