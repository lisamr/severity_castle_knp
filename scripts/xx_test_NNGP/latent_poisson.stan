functions{
  real nngp_w_lpdf(vector w, real sigmasq, real phi, matrix NN_dist,
  matrix NN_distM, int[,] NN_ind, int N, int M){
    vector[N] V;
    vector[N] I_Aw = w;
    int dim;
    int h;
    for (i in 2:N) {
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNdistM;
      matrix[ i < (M + 1)? (i - 1) : M, i < (M + 1)? (i - 1): M]
      iNNCholL;
      vector[ i < (M + 1)? (i - 1) : M] iNNcorr;
      vector[ i < (M + 1)? (i - 1) : M] v;
      row_vector[i < (M + 1)? (i - 1) : M] v2;
      dim = (i < (M + 1))? (i - 1) : M;
      if(dim == 1){iNNdistM[1, 1] = 1;}
      else{
        h = 0;
        for (j in 1:(dim - 1)){
          for (k in (j + 1):dim){
            h = h + 1;
            iNNdistM[j, k] = exp(- phi * NN_distM[(i - 1), h]);
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for(j in 1:dim){
          iNNdistM[j, j] = 1;
        }
      }
      iNNCholL = cholesky_decompose(iNNdistM);
      iNNcorr = to_vector(exp(- phi * NN_dist[(i - 1), 1:dim]));
      v = mdivide_left_tri_low(iNNCholL, iNNcorr);
      V[i] = 1 - dot_self(v);
      v2 = mdivide_right_tri_low(v', iNNCholL);
      I_Aw[i] = I_Aw[i] - v2 * w[NN_ind[(i - 1), 1:dim]];
    }
    V[1] = 1;
    return - 0.5 * ( 1 / sigmasq * dot_product(I_Aw, (I_Aw ./ V)) +
    sum(log(V)) + N * log(sigmasq));
  }
}
data {
  int<lower=1> N; // number of samples
  int<lower=1> M; // number of neighbors
  int<lower=1> P; // number of predictors
  vector[N] Y; // response
  matrix[N, P + 1] X; // predictor matrix + intercept
  int NN_ind[N - 1, M]; // j index matrix for M neighbors where j < i
  matrix[N - 1, M] NN_dist; // distances to those neighbors
  matrix[N - 1, (M * (M - 1)) %/% 2] NN_distM; //lower triangular part of distance matrix
  // define constants used in priors
  vector[P+1] uB; // mean vector for MVN beta prior
  cov_matrix[P+1] VB; // cov matrix for MVN beta prior
  real<lower=0> ss; // scale par for sigma
  real<lower=0> st; // scale par for tau
  real<lower=0> ap; // shape par for phi
  real<lower=0> bp; // scale par for phi
}
transformed data {
  cholesky_factor_cov[P + 1] L_VB;
  L_VB = cholesky_decompose(VB);
}
parameters{
  vector[P + 1] beta; // fixed effects
  real<lower = 0> sigma; // marginal SD for GP 
  real<lower = 0> tau; // noise for likelihood
  real<lower = 0> phi; // decay parameter for GP
  vector[N] w; 
}
transformed parameters {
  real sigmasq = sigma^2;
  real tausq = tau^2;
}
model{
  beta ~ multi_normal_cholesky(uB, L_VB);
  phi ~ gamma(ap, bp);
  sigma ~ normal(0, ss);
  tau ~ normal(0, st);
  w ~ nngp_w(sigmasq, phi, NN_dist, NN_distM, NN_ind, N, M);
  Y ~ normal(X * beta + w, tau);
}
