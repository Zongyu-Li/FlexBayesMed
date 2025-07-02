functions {
  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real tau, real sigma, real gamma) {
    int n = num_elements(y);
    real c = sqrt(pi() / 2) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] + sigma * c >= 0) {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i] + sigma * c) / (sigma * gamma) | 0, 1);
      } else {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i] + sigma * c) / (sigma * (1 / gamma)) | 0, 1);
      }
    }

    return log_l1 + sum(log_y_comp);
  }
}

data {
  int<lower=0> nY; // number of observations
  vector[nY] Y;    // observations
  vector[nY] M;    // observations
  vector[nY] X;    // predictor values

  // Parameters for gamma's prior
  real<lower=0> alpha_gamma;
  real<lower=0> beta_gamma;
  real<lower=0> gamma_lower; // left trunc
  real<lower=0> gamma_upper; // right trunc
}

parameters {
  real intercept; // intercept
  real tau;
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
}

model {
  // Prior: p(intercept, tau, sigma) ∝ 1/sigma
  target += -log(sigma);

  // Prior: gamma
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
          - log_diff_exp(
              gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
              gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
            );

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, tau, sigma, gamma);
}

