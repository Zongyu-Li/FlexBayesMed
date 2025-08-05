functions {
  // Computes m(γ,ν) = f(ν) * (γ - γ⁻¹) for mean adjustment
  real mean_adjustment(real v) {
    real log_num = log(2) + log(v) + lgamma((v + 1)/2);
    real log_denom = 0.5*(log(pi()) + log(v)) + log(v - 1) + lgamma(v/2);
    return exp(log_num - log_denom);
  }

  // MATPT log-likelihood for regression
  real matpt_lpdf(
    vector y,           // Response vector
    matrix X,           // Design matrix
    real intercept,
    vector beta,
    real sigma,
    real gamma,
    real v
  ) {
    int n = size(y);
    real gamma_inv = 1/gamma;
    real m_adj = mean_adjustment(v) * (gamma - gamma_inv);  // m(γ,ν)
    real log_const = n * (log(2) - log(sigma) - log(gamma + gamma_inv));
    vector[n] lp_terms;

    for (i in 1:n) {
      real mu = intercept + dot_product(X[i], beta);
      real resid = y[i] - mu;                   // Raw residual
      real shifted_resid = resid + sigma * m_adj;  // ε + m(γ,ν)

      if (shifted_resid >= 0) {
        // Right side: scale = σγ
        real scaled_resid = shifted_resid / (sigma * gamma);
        lp_terms[i] = student_t_lpdf(scaled_resid | v, 0, 1);
      } else {
        // Left side: scale = σ/γ
        real scaled_resid = shifted_resid / (sigma * gamma_inv);
        lp_terms[i] = student_t_lpdf(scaled_resid | v, 0, 1);
      }
    }
    return log_const + sum(lp_terms);
  }
}

data {
  // Sample size and predictors
  int<lower=0> n;       // Number of observations
  int<lower=1> k;       // Number of predictors

  // Response vector
  vector[n] y;
  // Predictor matrix (n x k, no intercept column)
  matrix[n, k] X;

  // MATPT hyperparameters
  real<lower=0> alpha_gamma;    // Gamma shape
  real<lower=0> beta_gamma;     // Gamma rate
  real<lower=0> gamma_lower;    // Gamma lower bound
  real<lower=0> gamma_upper;    // Gamma upper bound
  real<lower=0> d_v;            // Exponential rate for v
  real<lower=1> v_lower;        // v lower bound (enforces ν > 1)
}

parameters {
  // Regression parameters
  real intercept;
  vector[k] beta;
  real<lower=0> sigma;

  // MATPT shape parameters
  real<lower=gamma_lower, upper=gamma_upper> gamma;
  real<lower=v_lower> v;
}

model {
  //---------- Priors ----------//
  target += -log(sigma);

  // Gamma prior for skewness parameter (truncated)
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
            - log_diff_exp(
                gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
                gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
              );

  // Exponential prior for degrees of freedom (truncated)
  target += exponential_lpdf(v | d_v)
            - exponential_lccdf(v_lower | d_v);

  //---------- Likelihood ----------//
  target += matpt_lpdf(y | X, intercept, beta, sigma, gamma, v);
}
