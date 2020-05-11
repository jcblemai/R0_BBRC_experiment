// Bayesian mixture model for Accelerated Failure Time survival models with competing risk
// This version uses Gamma Distributions
data {
  int<lower=1> N;       // number of observations                                     
  int<lower=1> M;       // number of different events (in addition to censored)
  real<lower=0> times[N];          // times to events
 // int<lower=0, upper=M> events[N]; // label of events, 0 indicates censored;
  int<lower=0> indices[M+1, N];    // indices of events corresponding to each eveny type, first row is censored observations
  int<lower=0> N_events[M];        // number of events in each category, used to access matrix of indices
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
    target += gamma_lpdf(times[indices[m, 1:N_events[m]]] | alpha[m], alpha[m]/mu[m]) + log(lambda[m]); 
  }
  
  // Priors on mu and alpha
  mu ~ cauchy(0, 2.5);
  alpha ~ cauchy(0,2.5);   
}
