/* Response NNGP model */
functions{
  real nngp_w_lpdf(vector w_b1, real sigmasq, real phi, matrix NN_dist,
  matrix NN_distM, int[,] NN_ind, int N, int M,
  real intercept){
    
    vector[N] V;
    vector[N] w = w_b1 - intercept;
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
      iNNcorr = to_vector(exp(- phi * NN_dist[(i - 1), 1: dim]));
      
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
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> P;
  vector[N] Y;
  matrix[N, P + 1] X;
  int NN_ind[N - 1, M];
  matrix[N - 1, M] NN_dist;
  matrix[N - 1, (M * (M - 1)) %/% 2] NN_distM;
  real se;
  real ss;
  real ap;
  real bp;
}
parameters{
  vector[P + 1] beta;
  real<lower = 0> eta;
  real<lower = 0> sigma;
  real<lower = 0> rho;
  vector[N] w_b1;
}
transformed parameters {
  real etasq = square(eta);
  real sigmasq = square(sigma);
}
model{
  beta ~ normal(0, 1);
  rho ~  inv_gamma(ap, bp); //gamma(ap, bp);
  eta ~ normal(0, se);
  sigma ~ normal(0, ss);
  w_b1 ~ nngp_w(etasq, rho, NN_dist, NN_distM, NN_ind, N, M, beta[1]);
  Y ~ normal(block(X, 1, 2, N, P) * tail(beta, P) + w_b1, sigma);
}
