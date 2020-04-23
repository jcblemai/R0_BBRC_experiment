data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int<lower=1> C;          // number of integer shape values to evaluate
  real y[N];               // observations
  real<lower=0, upper=1> C_prior[C];
  real<lower=0> mu_scale[max(K, 2)];
}
parameters {
  simplex[K] theta;           // mixing proportions
  positive_ordered[K] mu;  // mean of mixture components
}
transformed parameters {
  matrix[K,C] lp [N] ;  // cache log calculation
  
  for (n in 1:N) {
    for (k in 1:K) {
      for (c in 1:C) {
        lp[n][k,c] = log(theta[k]) + log(C_prior[c]);
      }
    }
  }
  
  for (n in 1:N) {
    for (c in 1:C) {
      for (k in 1:K)
      lp[n][k, c] += gamma_lpdf(y[n] | c, c/mu[k]);
    }
  }
}
model {
  
  for(k in 1:K)
    mu[k] ~ cauchy(mu_scale[k], 2.5);
  
  for (n in 1:N) {
    target += log_sum_exp(to_vector(lp[n]));
  }
}
generated quantities {
  vector[K] ll_k[N];
  vector[C] ll_c[K];
  vector[N] log_lik;
  
  for (n in 1:N)
  log_lik[n] = 0.0;
  for (c in 1:C)
  for (k in 1:K)
  ll_c[k][c] = 0.0;
  
  for (n in 1:N)
  for (k in 1:K)
  ll_k[n][k] = 0.0;
  
  for (n in 1:N)
  log_lik[n] += log_sum_exp(to_vector(lp[n]));
  
  for (c in 1:C)
  for (k in 1:K)
  for (n in 1:N)
  ll_c[k][c] += lp[n][k, c];
  
  for (k in 1:K)
  ll_c[k] = softmax(ll_c[k]);
  
  for (n in 1:N)
  for (k in 1:K)
  ll_k[n][k] += log_sum_exp(lp[n, k,]);
  
  for (n in 1:N)
  ll_k[n] = softmax(ll_k[n]);
}

