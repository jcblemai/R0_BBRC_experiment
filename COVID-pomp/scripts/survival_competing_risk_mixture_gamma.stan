// Bayesian mixture model for Accelerated Failure Time survival models with competing risk
// Based on: http://proceedings.mlr.press/v84/bellot18a.html
// This version uses Gamma Distributions
data {
  int<lower=1> N;       // number of observations                                     
  int<lower=1> M;       // number of different events (in addition to censored)
  real<lower=0> times[N];          // times to events
  // int<lower=0, upper=M> events[N]; // label of events, 0 indicates censored;
  int<lower=0> N_events[M+1];        // number of events in each category, used to access matrix of indices
  int<lower=0> indices[M+1, max(N_events)];    // indices of events corresponding to each eveny type, first row is censored observations
  int<lower=1> K;          // number of mixture components
}
parameters {
  real<lower=0> alpha[M, K];  // shape parameter of gamma distributions                                     
  real<lower=0> mu[M, K];     // mean of gamma distributions (rate parameter = alpha/mu)  
  simplex[M] lambda;          // probability of different events
  simplex[K] theta;           // mixing proportions
}
transformed parameters {
  vector[K] lp [N] ;  // cache log calculation
  
  for (i in 1:N_events[1]) {
    for (k in 1:K) {
      lp[indices[1, i]][k] = 1;
      for (m in 1:M) {
        lp[indices[1, i]][k] -= gamma_cdf(times[indices[1, i]], alpha[m, k], alpha[m, k]/mu[m, k]) * lambda[m];
      }
      lp[indices[1, i]][k] = log(lp[indices[1, i]][k]) + log(theta[k]);
    }
  }
  
  // uncensored data 
  for (m in 1:M) {
    for (i in 1:N_events[m+1]) {
      for (k in 1:K) {
        lp[indices[m+1, i]][k] = log(lambda[m]) + gamma_lpdf(times[indices[m+1, i]] | alpha[m, k], alpha[m, k]/mu[m, k]) + log(theta[k]);
      }
    }
  }
}
model {
  for (n in 1:N) {
    target += log_sum_exp(to_vector(lp[n]));
  }
  
  // Priors on mu and alpha
  for (k in 1:K) {
    mu[, k] ~ cauchy(0, 2.5);
    alpha[, k] ~ cauchy(0,2.5);   
  }
}
generated quantities {
  vector[N] log_lik;
  
  for (n in 1:N)
  log_lik[n] = 0.0;
  
  for (n in 1:N)
  log_lik[n] += log_sum_exp(to_vector(lp[n]));
  
}
