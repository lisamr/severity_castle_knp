data{
  int N; // samples
  vector[N] dd_var; // mean of d-d var in each ring
  vector[N] y; // response variable
}
parameters{
  real a0;
  real beta;
  real<lower=0> sigma;
}
transformed parameters{
  vector[N] mu;
  
  
  // main model
  mu = a0 + beta * dd_var;
  
}
model{
  a0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  target += exponential_lpdf(sigma | 1);
  //sigma ~ exponential(1);
  
  y ~ normal(mu, sigma);
}
generated quantities{
  vector[N] y_rep;
  for(i in 1:N){
    y_rep[i] = normal_rng(mu[i], sigma);
  }
}
