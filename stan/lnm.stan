functions {
  vector phi_inv(row_vector mu, int K) {
    vector[K] e_mu;

    e_mu = exp(to_vector(mu));
    return append_row(e_mu / (1 + sum(e_mu)), 1 / (1 + sum(e_mu)));
  }

  vector sanitize_p(vector p, int K) {
    vector[K] p_new;
    int nan_ix;
    p_new = p;

    nan_ix = -1;
    for (k in 1:K) {
      if (is_nan(p_new[k])) {
        nan_ix = k;
      }
    }

    if (nan_ix != -1) {
      for (k in 1:K) {
        p_new[k] = .01 / (K - 1);
      }
      p_new[nan_ix] = .99;
    }

    return p_new;
  }
}
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] trt;
  array[N, K + 1] int<lower=0> y;
}

parameters {
  vector[K] beta1_t;
  vector[K] beta1_intercept;
  vector<lower=0>[K] sigmas1;
  matrix[N, K] mu;
}

model {
  // prior
  for (k in 1:K) {
    beta1_intercept[k] ~ normal(0, 1);
    sigmas1[k] ~ inv_gamma(50, 25);
    beta1_t[k] ~ normal(0, 1);
    
  }

  // likelihood
  vector[K + 1] xi;
  vector[K] e_mu;
  for (i in 1:N) {
    for (k in 1:K) {
      mu[i][k] ~ normal(beta1_intercept[k] + beta1_t[k] * trt[i], sigmas1[k]);

    }
    e_mu = to_vector(exp(mu[i]));
    xi = append_row(e_mu / (1 + sum(e_mu)), 1 / (1 + sum(e_mu)));
    y[i] ~ multinomial(to_vector(xi));
  }
}

generated quantities {
  // simulate samples from a full posterior
  vector[K + 1] p_sim;
  matrix[N, K] mu_sim;
  array[N, K + 1] int<lower=0> y_sim;
  array[N, K + 1] real<lower=0> p_sim_arr;

  for (i in 1:N) {
    for (k in 1:K) {
      mu_sim[i, k] = normal_rng(beta1_intercept[k] + beta1_t[k] * trt[i], sigmas1[k]);
    }
    p_sim = phi_inv(mu_sim[i], K);
    p_sim = sanitize_p(p_sim, K + 1);
    y_sim[i] = multinomial_rng(to_vector(p_sim), sum(y[i]));

    // obtain simulated relative abundance
    for (k in 1:(K + 1)) {
      p_sim_arr[i][k] = y_sim[i][k] * 1.0 / sum(y_sim[i]);
    }
  }

  // counterfactuals
  array[2, N, K + 1] real<lower=0> p_c; // conterfactual in probability scale
  array[2, N, K + 1] int<lower=0> y_c;
  vector[K + 1] effect;
  matrix[N, K] mu_c;
  vector[K + 1] temp_pc;
  matrix[2, K + 1] sum_p; 
  for (t in 0:1) {
    for (k in 1:(K + 1)){
      sum_p[t + 1][k] = 0;
    }  
    for (i in 1:N) {
      for (k in 1:K) {
        mu_c[i, k] = normal_rng(beta1_intercept[k] + beta1_t[k] * t, sigmas1[k]);
        
        temp_pc = phi_inv(mu_c[i], K);
        temp_pc = sanitize_p(temp_pc, K + 1);
        y_c[t + 1][i] = multinomial_rng(to_vector(temp_pc), sum(y[i]));
      }// loop through k

      for (k in 1:(K + 1)){
          
          sum_p[t + 1][k] += temp_pc[k];

          // counterfactuals in probability scale
          p_c[t + 1][i][k] = y_c[t + 1][i][k] * 1.0 / sum(y[i]);
      }  
    }// loop through i
  }// loop through t
  for (k in 1:(K + 1)){
    effect[k] = (sum_p[2, k] - sum_p[1, k]) / 1.0 / N;
  }  
  
}



