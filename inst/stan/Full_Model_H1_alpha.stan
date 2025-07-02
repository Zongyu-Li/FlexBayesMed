functions {
  real f(real nu) {
    real log_numerator = log(2) + log(nu) + lgamma((nu + 1) / 2);
    real log_denominator = 0.5 * (log(pi()) + log(nu)) + log(nu - 1) + lgamma(nu / 2);
    real log_ratio = log_numerator - log_denominator;
    return exp(log_ratio);
  }

  real log_likelihood_lpdf(vector m, vector x, real intercept, real alpha, real sigma, real gamma, real nu) {
    int n = num_elements(m);
    real c = f(nu) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept - alpha * x[i] + sigma * c >= 0) {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept - alpha * x[i] + sigma * c) / (sigma * gamma) | nu, 0, 1);
      } else {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept - alpha * x[i] + sigma * c) / (sigma * (1 / gamma)) | nu, 0, 1);
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

  // Parameters for nu's prior
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept; // intercept
  real alpha; // slope
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, alpha, sigma) ∝ 1/sigma
  target += -log(sigma);

  // Prior: gamma
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
          - log_diff_exp(
              gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
              gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
            );

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, alpha, sigma, gamma, nu);
}

