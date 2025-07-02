functions {
  real log_likelihood_lpdf(vector m, vector x, real intercept, real sigma) {
    int n = num_elements(m);
    real log_l1 = - n * log(sigma);
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept>= 0) {
        log_m_comp[i] = normal_lpdf((m[i] - intercept) / (sigma) | 0, 1);
      } else {
        log_m_comp[i] = normal_lpdf((m[i] - intercept) / (sigma) | 0, 1);
      }
    }
    
    return log_l1 + sum(log_m_comp);
  }
}

data {
  int<lower=0> nM; // number of observations
  vector[nM] M;    // observations
  vector[nM] X;    // predictor values
}

parameters {
  real intercept; // intercept
  real<lower=0> sigma; // standard deviation
}

model {
  // Prior: p(intercept, sigma) ∝ 1/sigma
  target += -log(sigma);
  
  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, sigma);
}

