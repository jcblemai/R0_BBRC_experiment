data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
  real<lower=0> mu_scale[max(K, 2)];
}
parameters {
  simplex[K] theta;           // mixing proportions
  positive_ordered[K] mu;  // mean of mixture components
  real<lower=0, upper=10> alpha [K];
}
transformed parameters {
  vector[K] lp [N] ;  // cache log calculation
  
  for (n in 1:N) {
    for (k in 1:K) {
        lp[n][k] = log(theta[k]);
    }
  }
  
  for (n in 1:N) {
      for (k in 1:K)
      lp[n][k] += gamma_lpdf(y[n] | alpha[k], alpha[k]/mu[k]);
  }
}
model {
  for(k in 1:K)
    mu[k] ~ cauchy(mu_scale[k], 2.5);
    
  alpha ~ cauchy(0,2.5);
  for (n in 1:N) {
    target += log_sum_exp(to_vector(lp[n]));
  }
}
generated quantities {
  vector[K] ll_k[N];
  vector[N] log_lik;
  
  for (n in 1:N)
  log_lik[n] = 0.0;
  
  for (n in 1:N)
  for (k in 1:K)
  ll_k[n][k] = 0.0;
  
  for (n in 1:N)
  log_lik[n] += log_sum_exp(to_vector(lp[n]));
  
  for (n in 1:N)
  ll_k[n] = softmax(lp[n,]);
}

