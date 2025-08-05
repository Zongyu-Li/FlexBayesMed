functions {
  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real beta, real tau, real sigma) {
    int n = num_elements(y);
    real log_l1 = - n * log(sigma);
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] - beta * m[i]>= 0) {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i] - beta * m[i]) / (sigma) | 0, 1);
      } else {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i] - beta * m[i]) / (sigma) | 0, 1);
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
}

parameters {
  real intercept; // intercept
  real beta; // slope
  real tau;
  real<lower=0> sigma; // standard deviation
}

model {
  // Prior: p(intercept, beta, tau, sigma) ∝ 1/sigma
  target += -log(sigma);

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, beta, tau, sigma);
}

