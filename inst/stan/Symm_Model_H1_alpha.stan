functions {
  real log_likelihood_lpdf(vector m, vector x, real intercept, real alpha, real sigma, real nu) {
    int n = num_elements(m);
    real log_l1 = - n * log(sigma);
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept - alpha * x[i] >= 0) {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept - alpha * x[i]) / (sigma) | nu, 0, 1);
      } else {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept - alpha * x[i]) / (sigma) | nu, 0, 1);
      }
    }

    return log_l1 + sum(log_m_comp);
  }
}

data {
  int<lower=0> nM; // number of observations
  vector[nM] M;    // observations
  vector[nM] X;    // predictor values

  // Parameters for nu's prior
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept; // intercept
  real alpha; // slope
  real<lower=0> sigma; // standard deviation
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, alpha, sigma) ∝ 1/sigma
  target += -log(sigma);

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, alpha, sigma, nu);
}

