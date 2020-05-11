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
}
parameters {
  real<lower=0> alpha[M];    // shape parameter of gamma distributions                                     
  real<lower=0> mu[M];       // mean of gamma distributions (rate parameter = alpha/mu)  
  simplex[M] lambda;          // probability of different events
}
model {
  // censored data
  // cash computation of censored observations survival functions
  real Si[N_events[1]];
  
  for (i in 1:N_events[1]) {
    Si[i] = 1;
    for (m in 1:M) {
      Si[i] -= gamma_cdf(times[indices[1, i]], alpha[m], alpha[m]/mu[m]) * lambda[m];
    }
    target += log(Si[i]);
  }
  // uncensored data 
  for (m in 1:M) {
    target += gamma_lpdf(times[indices[m+1, 1:N_events[m+1]]] | alpha[m], alpha[m]/mu[m]);
    target += N_events[m+1] * log(lambda[m]); 
  }
  
  // Priors on mu and alpha
  mu ~ cauchy(0, 2.5);
  alpha ~ cauchy(0,2.5);   
}
generated quantities {
  vector[N] log_lik; 
  
  for (i in 1:N_events[1]) {
    log_lik[indices[1, i]] = 1;
    for (m in 1:M) {
      log_lik[indices[1, i]] -= gamma_cdf(times[indices[1, i]], alpha[m], alpha[m]/mu[m]) * lambda[m];
    }
    log_lik[indices[1, i]] = log(log_lik[indices[1, i]]);
  }
  // uncensored data 
  for (m in 1:M) {
    for (i in 1:N_events[m+1]) {
      log_lik[indices[m+1, i]] = log(lambda[m]) + gamma_lpdf(times[indices[m+1, i]] | alpha[m], alpha[m]/mu[m]);
    }
  }
}
