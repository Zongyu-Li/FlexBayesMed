functions {
  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real beta, real tau, real sigma, real nu) {
    int n = num_elements(y);
    real log_l1 = - n * log(sigma);
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] - beta * m[i] >= 0) {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i] - beta * m[i]) / (sigma) | nu, 0, 1);
      } else {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i] - beta * m[i]) / (sigma) | nu, 0, 1);
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

  // Parameters for nu's prior
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept; // intercept
  real beta; // slope
  real tau;
  real<lower=0> sigma; // standard deviation
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, beta, tau, sigma) ∝ 1/sigma
  target += -log(sigma);

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, beta, tau, sigma, nu);
}

