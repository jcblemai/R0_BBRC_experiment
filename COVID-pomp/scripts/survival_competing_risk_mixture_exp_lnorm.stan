// Bayesian mixture model for Accelerated Failure Time survival models with competing risk
// Based on: http://proceedings.mlr.press/v84/bellot18a.html
// This version uses a mixture of an exponential and a Log-Normal Distribution
data {
  int<lower=1> N;       // number of observations                                     
  int<lower=1> M;       // number of different events (in addition to censored)
  real<lower=0> times[N];          // times to events
  // int<lower=0, upper=M> events[N]; // label of events, 0 indicates censored;
  int<lower=0> N_events[M+1];        // number of events in each category, used to access matrix of indices
  int<lower=0> indices[M+1, max(N_events)];    // indices of events corresponding to each eveny type, first row is censored observations
  int<lower=1> K[M];
}
parameters {
  real<lower=0> sigma[M];  //                                    
  real<lower=0> mu[M];     // 
  real<lower=0, upper=1> mu_diff[M];  // to ensure mu_exp < mu
  simplex[M] lambda;          // probability of different events
  real<lower=0, upper=1> theta[M];           // mixing proportions
}
transformed parameters {
  vector[max(K)] lp [N] ;  // cache log calculation
  real mu_exp[M];
  real thetak;
  thetak = 1;
  for (m in 1:M)
  mu_exp[m] = exp(mu[m] + (sigma[m]^2)/2) * mu_diff[m]; // as fraction of lognorm mean
  
  for (i in 1:N_events[1]) {
    for (k in 1:max(K)) {
      lp[indices[1, i]][k] = 0;
    }
  }
  
  // Censored data
  for (i in 1:N_events[1]) {
    for (m in 1:M) {
      for (k in 1:K[m]) {
        if(K[m] > 1){
          if (k == 1) {
            thetak =  theta[m];
          } else {
            thetak = 1-theta[m];
          }
          if (k == 1) {
            lp[indices[1, i]][k] += exponential_cdf(times[indices[1, i]], 1/mu_exp[m]) * lambda[m] * thetak;
          } else {
            lp[indices[1, i]][k] += lognormal_cdf(times[indices[1, i]], mu[m], sigma[m]) * lambda[m] * thetak;
          }
        } else {
          lp[indices[1, i]][k] += lognormal_cdf(times[indices[1, i]], mu[m], sigma[m]) * lambda[m];
        }
      }
    }
  }
  
  // uncensored data 
  for (m in 1:M) {
    for (i in 1:N_events[m+1]) {
      for (k in 1:K[m]) {
        if(K[m] > 1){
          if (k == 1) {
            thetak =  theta[m];
          } else {
            thetak = 1-theta[m];
          }
          if (k == 1) {
            lp[indices[m+1, i]][k] = log(lambda[m]) + exponential_lpdf(times[indices[m+1, i]] | 1/mu_exp[m]) + log(thetak);
          } else  {
            lp[indices[m+1, i]][k] = log(lambda[m]) + lognormal_lpdf(times[indices[m+1, i]] | mu[m], sigma[m]) + log(thetak);
          }
        } else {
          lp[indices[m+1, i]][k] = log(lambda[m]) + lognormal_lpdf(times[indices[m+1, i]] | mu[m], sigma[m]);
          lp[indices[m+1, i]][k+1] = -1e6; // numerical trick to set prob to 0 for missing group
        }
      }
    }
  }
}
model {
  // Censored observations
  for (i in 1:N_events[1]) {
    target += log(1-sum(to_vector(lp[indices[1, i]])));
  }
  
  // Uncensored
  for (m in 1:M) {
    for (i in 1:N_events[m+1]) {
      target += log_sum_exp(to_vector(lp[indices[m+1, i]]));
    }
  }

  // Priors on mu and alpha
  mu_diff ~ std_normal();
  mu ~ cauchy(0, 1);
  sigma ~ cauchy(0, 1);   
}
generated quantities {
  vector[N] log_lik;
  
  // Censored observations
  for (i in 1:N_events[1]) {
    log_lik[indices[1, i]] = log(1-sum(to_vector(lp[indices[1, i]])));
  }
  
  // Uncensored
  for (m in 1:M) {
    for (i in 1:N_events[m+1]) {
      log_lik[indices[m+1, i]] = log_sum_exp(to_vector(lp[indices[m+1, i]]));
    }
  }
}
