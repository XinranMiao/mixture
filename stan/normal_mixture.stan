data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] trt;
  matrix[N, K] z;
}

parameters {
  vector[K] beta1_t;
  vector[K] beta1_intercept;
  vector[K] beta0_intercept;
  vector<lower=0>[K] sigmas1;
  vector<lower=0>[K] sigmas0;
  vector<lower=0,upper=1>[K] p1;
}

model {
  // prior
  for (k in 1:K) {
    beta1_intercept[k] ~ normal(0, 1);
    beta0_intercept[k] ~ normal(0, 1);
    sigmas1[k] ~ inv_gamma(50, 25);
    sigmas0[k] ~ inv_gamma(50, 25);
    beta1_t[k] ~ normal(0, 1);
    
  }

  // likelihood
  for (i in 1:N) {
    for (k in 1:K) {
      target += log_mix( p1[k],
                    normal_lpdf(z[i][k] | beta1_intercept[k] + beta1_t[k] * trt[i], sigmas1[k]),
                    normal_lpdf(z[i][k] | beta0_intercept[k], sigmas0[k])) ;

    }
  }
}


generated quantities {
  // simulate samples from a full posterior
  matrix[N, K] z_sim;
  for (i in 1:N) {
    for (k in 1:K) {
        int delta = bernoulli_rng(p1[k]);
        if (delta == 1) {
          z_sim[i, k] = normal_rng(beta1_intercept[k] + beta1_t[k] * trt[i], sigmas1[k]);
        } else {
          z_sim[i, k] = normal_rng(beta0_intercept[k], sigmas0[k]);
        }
      
    }
  }
}
