functions {
  real log_likelihood_lpdf(vector m, vector x, real intercept, real sigma, real gamma) {
    int n = num_elements(m);
    real c = sqrt(pi() / 2) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept + sigma * c >= 0) {
        log_m_comp[i] = normal_lpdf((m[i] - intercept + sigma * c) / (sigma * gamma) | 0, 1);
      } else {
        log_m_comp[i] = normal_lpdf((m[i] - intercept + sigma * c) / (sigma * (1 / gamma)) | 0, 1);
      }
    }

    return log_l1 + sum(log_m_comp);
  }
}

data {
  int<lower=0> nM; // number of observations
  vector[nM] M;    // observations
  vector[nM] X;    // predictor values

  // Parameters for gamma's prior
  real<lower=0> alpha_gamma;
  real<lower=0> beta_gamma;
  real<lower=0> gamma_lower; // left trunc
  real<lower=0> gamma_upper; // right trunc
}

parameters {
  real intercept; // intercept
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
}

model {
  // Prior: p(intercept, sigma) ∝ 1/sigma
  target += -log(sigma);

  // Prior: gamma
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
          - log_diff_exp(
              gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
              gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
            );

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, sigma, gamma);
}

