functions {
  vector stan_interp(real SS, matrix E) {
    int index = cols(E);
    
    // If SS is already above the maximum index then
    // we can skip the interpolation
    if (index < SS) {
      return E[:, index];
    }
    
    while (index > SS) {
      index -= 1;
    }
    
    real prop_scales = SS - index;
    
    // Assumes same value of SS for every row
    return E[:,index] * (1 - prop_scales) + E[:,index + 1] * prop_scales;
  }
}
data {
  int<lower=0> N; // number of sites 
  int<lower=1> M; // number of buffers
  vector[N] y; // observations
  matrix[N, M] dd_var; // covariate values within each buffer
}
transformed data{
  real off = 0.00001;
}
parameters {
  real a0;
  real beta;
  real<lower=0> sigma; // standard deviation
  real<lower=1, upper=M-off> ss; // scale selection parameter
}
transformed parameters {
  vector[N] w = stan_interp(ss, dd_var);
}
model {
  // priors
  a0 ~ normal(0, 2);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 2);
  ss ~ lognormal(1, 1) T[,M-off]; 
  
  y ~ normal(a0 + beta * w, sigma);
  
}
