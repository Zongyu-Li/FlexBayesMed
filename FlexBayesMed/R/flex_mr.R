#' Fit a Linear Regression with Centred Two-Piece Student \eqn{t} (CTPT) Errors
#'
#' @param X A numeric matrix of predictors (\eqn{n \times k}).
#' @param Y A numeric vector of length \eqn{n} (the response).
#' @param model A string: which model to use ("full", "gamma-only", "nu-only", "normal")
#' @param alpha_gamma,beta_gamma Shape and rate for the Gamma prior on skewness parameter \eqn{\gamma}.
#' @param gamma_lower,gamma_upper Truncation bounds for \eqn{\gamma} (\eqn{0 < \gamma_{\text{lower}} \leq \gamma \leq \gamma_{\text{upper}}}).
#' @param d_nu Rate for the exponential prior on degrees-of-freedom \eqn{\nu}.
#' @param nu_lower Lower truncation bound for \eqn{\nu} (must satisfy \eqn{\nu > 1}).
#' @param iter,warmup,chains,seed,verbose,refresh,... Additional arguments passed to `rstan::stan()`.
#'
#' @return An object of class `stanfit`.
#' @importFrom rstan stan
#' @export
#'
flex_mr <- function(X, Y, model = "full", alpha_gamma = 2, beta_gamma = 2, gamma_lower = 0.05,
                    gamma_upper = 20, d_nu = 0.01, nu_lower = 1.01, iter = 30000, warmup = floor(iter/2),
                    chains = 1, seed = 123, verbose = TRUE, refresh = 2500, ...) {

  valid_models <- c("full", "gamma-only", "nu-only", "normal")
  if (!model %in% valid_models) {
    stop("Argument 'model' must be one of: 'full', 'gamma-only', 'nu-only', 'normal'.")
  }

  if(model=="full"){

  ## Construct data list for Stan
  stan_data <- list(
    n = nrow(X),
    k = ncol(X),
    y = as.vector(Y),
    X = X,
    alpha_gamma = alpha_gamma,
    beta_gamma  = beta_gamma,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    d_nu = d_nu,
    nu_lower = nu_lower
  )

  ## Stan model code (CTPT linear regression)
  model_code <- "
functions {
  // Computes m(gamma,nu) = f(nu) * (gamma - 1/gamma) for mean adjustment
  real mean_adjustment(real nu) {
    real log_num = log(2) + log(nu) + lgamma((nu + 1)/2);
    real log_denom = 0.5*(log(pi()) + log(nu)) + log(nu - 1) + lgamma(nu/2);
    return exp(log_num - log_denom);
  }

  // CTPT log-likelihood for regression
  real ctpt_lpdf(
    vector y,
    matrix X,
    real intercept,
    vector beta,
    real sigma,
    real gamma,
    real nu
  ) {
    int n = size(y);
    real gamma_inv = 1/gamma;
    real m_adj = mean_adjustment(nu) * (gamma - gamma_inv);  // m(gamma,nu)
    real log_const = n * (log(2) - log(sigma) - log(gamma + gamma_inv));
    vector[n] lp_terms;

    for (i in 1:n) {
      real mu = intercept + dot_product(X[i], beta);
      real resid = y[i] - mu;                  // Raw residual
      real shifted_resid = resid + sigma*m_adj;// resid + m(gamma,nu)

      if (shifted_resid >= 0) {
        // Right side: scale = sigma*gamma
        real scaled_resid = shifted_resid / (sigma * gamma);
        lp_terms[i] = student_t_lpdf(scaled_resid | nu, 0, 1);
      } else {
        // Left side: scale = sigma/gamma
        real scaled_resid = shifted_resid / (sigma / gamma);
        lp_terms[i] = student_t_lpdf(scaled_resid | nu, 0, 1);
      }
    }
    return log_const + sum(lp_terms);
  }
}

data {
  int<lower=0> n;
  int<lower=1> k;
  vector[n] y;
  matrix[n,k] X;

  // ctpt hyperparameters
  real<lower=0> alpha_gamma;
  real<lower=0> beta_gamma;
  real<lower=0> gamma_lower;
  real<lower=0> gamma_upper;
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept;
  vector[k] beta;
  real<lower=0> sigma;
  real<lower=gamma_lower, upper=gamma_upper> gamma;
  real<lower=nu_lower> nu;
}

model {
  // Jeffreys prior on sigma: p(sigma) 1/sigma
  target += -log(sigma);

  // Gamma prior for skewness (truncated)
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
            - log_diff_exp(
                gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
                gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
              );

  // Exponential prior for df (truncated)
  target += exponential_lpdf(nu | d_nu)
            - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += ctpt_lpdf(y | X, intercept, beta, sigma, gamma, nu);
}
"
  }

  if(model=="nu-only"){
    ## Construct data list for Stan
    stan_data <- list(
      n = nrow(X),
      k = ncol(X),
      y = as.vector(Y),
      X = X,
      d_nu = d_nu,
      nu_lower = nu_lower
    )

    model_code <- "
functions {
  // ctpt log-likelihood for regression
  real ctpt_lpdf(vector y, matrix X, real intercept, vector beta, real sigma, real nu) {
    int n = size(y);
    real log_const = n * (-log(sigma));
    vector[n] lp_terms;

    for (i in 1:n) {
      real mu = intercept + dot_product(X[i], beta);
      real resid = y[i] - mu;                  // Raw residual
      real shifted_resid = resid;
      real scaled_resid = shifted_resid / (sigma);
      lp_terms[i] = student_t_lpdf(scaled_resid | nu, 0, 1);

    }
    return log_const + sum(lp_terms);
  }
}

data {
  int<lower=0> n;
  int<lower=1> k;
  vector[n] y;
  matrix[n,k] X;

  // ctpt hyperparameters
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept;
  vector[k] beta;
  real<lower=0> sigma;
  real<lower=nu_lower> nu;
}

model {
  // Jeffreys prior on sigma: p(sigma) 1/sigma
  target += -log(sigma);

  // Exponential prior for df (truncated)
  target += exponential_lpdf(nu | d_nu)
            - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += ctpt_lpdf(y | X, intercept, beta, sigma, nu);
}
"

  }

  if(model=="gamma-only"){
    stan_data <- list(
      n = nrow(X),
      k = ncol(X),
      y = as.vector(Y),
      X = X,
      alpha_gamma = alpha_gamma,
      beta_gamma  = beta_gamma,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper
    )

    model_code <- "
functions {
  real mean_adjustment_gamma_only(real gamma) {
    return sqrt(2 / pi()) * (gamma - 1/gamma);
  }

  real twopiece_normal_lpdf(
    vector y,
    matrix X,
    real intercept,
    vector beta,
    real sigma,
    real gamma
  ) {
    int n = size(y);
    real log_const = n * (log(2) - log(sigma) - log(gamma + 1/gamma));
    vector[n] lp_terms;
    real c = mean_adjustment_gamma_only(gamma);

    for (i in 1:n) {
      real mu = intercept + dot_product(X[i], beta);
      real resid = y[i] - mu;
      real shifted_resid = resid + sigma * c;

      if (shifted_resid >= 0) {
        real scaled_resid = shifted_resid / (sigma * gamma);
        lp_terms[i] = normal_lpdf(scaled_resid | 0, 1) - log(gamma);
      } else {
        real scaled_resid = shifted_resid / (sigma / gamma);
        lp_terms[i] = normal_lpdf(scaled_resid | 0, 1) + log(gamma);
      }
    }
    return log_const + sum(lp_terms);
  }
}

data {
  int<lower=0> n;
  int<lower=1> k;
  vector[n] y;
  matrix[n,k] X;
  real<lower=0> alpha_gamma;
  real<lower=0> beta_gamma;
  real<lower=0> gamma_lower;
  real<lower=0> gamma_upper;
}

parameters {
  real intercept;
  vector[k] beta;
  real<lower=0> sigma;
  real<lower=gamma_lower, upper=gamma_upper> gamma;
}

model {
  target += -log(sigma);

  // Truncated gamma prior
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
            - log_diff_exp(
                gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
                gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
              );

  // likelihood
  target += twopiece_normal_lpdf(y | X, intercept, beta, sigma, gamma);
}
    "

  }

  if(model=="normal"){
    stan_data <- list(
      n = nrow(X),
      k = ncol(X),
      y = as.vector(Y),
      X = X,
      alpha_gamma = alpha_gamma,
      beta_gamma  = beta_gamma,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper
    )

    model_code <- "
functions {
  real my_normal_lpdf(
    vector y,
    matrix X,
    real intercept,
    vector beta,
    real sigma
  ) {
    int n = size(y);
    real log_const = n * (-log(sigma));
    vector[n] lp_terms;

    for (i in 1:n) {
      real mu = intercept + dot_product(X[i], beta);
      real resid = y[i] - mu;
      real scaled_resid = resid / sigma;
      lp_terms[i] = normal_lpdf(scaled_resid | 0, 1);
    }
    return log_const + sum(lp_terms);
  }
}

data {
  int<lower=0> n;
  int<lower=1> k;
  vector[n] y;
  matrix[n, k] X;
}

parameters {
  real intercept;
  vector[k] beta;
  real<lower=0> sigma;
}

model {
  // prior for scale: 1/sigma
  target += -log(sigma);

  // Normal likelihood
  target += my_normal_lpdf(y | X, intercept, beta, sigma);
}

    "
  }

## Run Stan
fit <- rstan::stan(
  model_code = model_code,
  data       = stan_data,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

return(fit)
}
