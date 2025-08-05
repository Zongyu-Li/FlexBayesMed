#' Mediation Analysis with Centred Two-Piece Student \eqn{t} (CTPT) Errors
#'
#' @param X A numeric vector of length \eqn{n} (the independent variable).
#' @param M A numeric vector of length \eqn{n} (the mediator).
#' @param Y A numeric vector of length \eqn{n} (the response).
#' @param model A string: which model to use ("full", "gamma-only", "nu-only", "normal")
#' @param q00 A non-negative numeric value giving the conditional probability
#'   \eqn{p(\alpha=0,\beta=0 \mid \alpha\beta=0)}. Defaults to \code{1/3}.
#' @param q01 A non-negative numeric value giving the conditional probability
#'   \eqn{p(\alpha=0,\beta\neq 0 \mid \alpha\beta=0)}. Defaults to \code{1/3}.
#'
#'   Here, \eqn{q00 + q01 + q10 = 1}, where \eqn{q10 = p(\alpha\neq 0,\beta=0 \mid \alpha\beta=0)}
#'   is computed internally as \code{1 - q00 - q01}. These three probabilities correspond
#'   to the disjoint ways in which the mediation effect can be zero (i.e.,
#'   \eqn{\alpha\beta = 0}). They must be non-negative and sum to 1.
#' @param alpha_gamma,beta_gamma Shape and rate for the Gamma prior on skewness parameter \eqn{\gamma}.
#' @param gamma_lower,gamma_upper Truncation bounds for \eqn{\gamma} (\eqn{0 < \gamma_{\text{lower}} \leq \gamma \leq \gamma_{\text{upper}}}).
#' @param d_nu Rate for the exponential prior on degrees-of-freedom \eqn{\nu}.
#' @param nu_lower Lower truncation bound for \eqn{\nu} (must satisfy \eqn{\nu > 1}).
#' @param iter,warmup,chains,seed,verbose,refresh,... Additional arguments passed to `rstan::stan()`.
#' @param repetitions,method,bridgesampling_cores,use_neff,maxiter,silent Additional arguments passed to `bridgesampling::bridge_sampler()`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{X2M_with_alpha}}{A \code{stanfit} object for the mediator model with \eqn{\alpha \neq 0}.}
#'   \item{\code{XM2Y_with_beta}}{A \code{stanfit} object for the outcome model with \eqn{\beta \neq 0}.}
#'   \item{\code{X2M_without_alpha}}{A \code{stanfit} object for the mediator model with \eqn{\alpha = 0}.}
#'   \item{\code{XM2Y_without_beta}}{A \code{stanfit} object for the outcome model with \eqn{\beta = 0}.}
#'   \item{\code{samples_med}}{A numeric vector of posterior draws for \eqn{\alpha \beta}.}
#'   \item{\code{bf_med}}{A numeric value giving the Bayes factor for the mediation effect.}
#' }
#' @importFrom rstan stan extract
#' @importFrom bridgesampling bridge_sampler bf
#' @export

flex_med <- function(X, M, Y, model = "full", q00 = 1/3, q01 = 1/3,
                     alpha_gamma = 2, beta_gamma = 2, gamma_lower = 0.05,
                    gamma_upper = 20, d_nu = 0.01, nu_lower = 1.01, iter = 30000,
                    warmup = floor(iter/2), chains = 1, seed = 123, verbose = FALSE,
                    refresh = 2500, repetitions = 1, method = "normal",
                    bridgesampling_cores = 1, use_neff = TRUE, maxiter = 1000, silent = FALSE,...) {

  if(q00 < 0 | q01 < 0){
    stop("The probabilities q00 and q01 should be both non-negative!")
  }
  q10 <- 1 - q00 - q01
  if(q10 < 0){
    stop("The sum of q00 and q01 should not exceed 1!")
  }

  valid_models <- c("full", "gamma-only", "nu-only", "normal")
  if (!model %in% valid_models) {
    stop("Argument 'model' must be one of: 'full', 'gamma-only', 'nu-only', 'normal'.")
  }


  ################## Full Model ########################
  if(model == "full"){
    model_code_H1_alpha <- "
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
  // Prior: p(intercept, alpha, sigma)  1/sigma
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
    "

    model_code_H1_beta <-"
functions {
  real f(real nu) {
    real log_numerator = log(2) + log(nu) + lgamma((nu + 1) / 2);
    real log_denominator = 0.5 * (log(pi()) + log(nu)) + log(nu - 1) + lgamma(nu / 2);
    real log_ratio = log_numerator - log_denominator;
    return exp(log_ratio);
  }

  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real beta, real tau, real sigma, real gamma, real nu) {
    int n = num_elements(y);
    real c = f(nu) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] - beta * m[i] + sigma * c >= 0) {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i] - beta * m[i] + sigma * c) / (sigma * gamma) | nu, 0, 1);
      } else {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i] - beta * m[i] + sigma * c) / (sigma * (1 / gamma)) | nu, 0, 1);
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

  // Parameters for nu's prior
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept; // intercept
  real beta; // slope
  real tau;
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, beta, tau, sigma)  1/sigma
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
  target += log_likelihood_lpdf(Y | M, X, intercept, beta, tau, sigma, gamma, nu);
}
    "

    model_code_H0_alpha <- "
functions {
  real f(real nu) {
    real log_numerator = log(2) + log(nu) + lgamma((nu + 1) / 2);
    real log_denominator = 0.5 * (log(pi()) + log(nu)) + log(nu - 1) + lgamma(nu / 2);
    real log_ratio = log_numerator - log_denominator;
    return exp(log_ratio);
  }

  real log_likelihood_lpdf(vector m, vector x, real intercept, real sigma, real gamma, real nu) {
    int n = num_elements(m);
    real c = f(nu) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept + sigma * c >= 0) {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept + sigma * c) / (sigma * gamma) | nu, 0, 1);
      } else {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept + sigma * c) / (sigma * (1 / gamma)) | nu, 0, 1);
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
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, sigma)  1/sigma
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
  target += log_likelihood_lpdf(M | X, intercept, sigma, gamma, nu);
}
    "

    model_code_H0_beta <- "
functions {
  real f(real nu) {
    real log_numerator = log(2) + log(nu) + lgamma((nu + 1) / 2);
    real log_denominator = 0.5 * (log(pi()) + log(nu)) + log(nu - 1) + lgamma(nu / 2);
    real log_ratio = log_numerator - log_denominator;
    return exp(log_ratio);
  }

  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real tau, real sigma, real gamma, real nu) {
    int n = num_elements(y);
    real c = f(nu) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] + sigma * c >= 0) {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i] + sigma * c) / (sigma * gamma) | nu, 0, 1);
      } else {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i] + sigma * c) / (sigma * (1 / gamma)) | nu, 0, 1);
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

  // Parameters for nu's prior
  real<lower=0> d_nu;
  real<lower=1> nu_lower;
}

parameters {
  real intercept; // intercept
  real tau;
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, tau, sigma)  1/sigma
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
  target += log_likelihood_lpdf(Y | M, X, intercept, tau, sigma, gamma, nu);
}
    "
    ## H1: Mediation exists
    print("Loading the Model: X2M_with_alpha")
    stan_data_H1_alpha <- list(
      nM = length(M),
      M = M,
      X = X,
      alpha_gamma = alpha_gamma,
      beta_gamma  = beta_gamma,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      d_nu = d_nu,
      nu_lower = nu_lower
    )

    fit_H1_alpha <- rstan::stan(
      model_code = model_code_H1_alpha,
      data       = stan_data_H1_alpha,
      iter       = iter,
      warmup     = warmup,
      chains     = chains,
      seed       = seed,
      verbose    = verbose,
      refresh    = refresh,
      ...
    )

    print("Loading the Model: XM2Y_with_beta")
    stan_data_H1_beta <- list(
      nY = length(Y),
      Y = Y,
      M = M,
      X = X,
      alpha_gamma = alpha_gamma,
      beta_gamma  = beta_gamma,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      d_nu = d_nu,
      nu_lower = nu_lower
    )

    fit_H1_beta <- rstan::stan(
      model_code = model_code_H1_beta,
      data       = stan_data_H1_beta,
      iter       = iter,
      warmup     = warmup,
      chains     = chains,
      seed       = seed,
      verbose    = verbose,
      refresh    = refresh,
      ...
    )

    ## H0: Mediation does not exist
    print("Loading the Model: X2M_without_alpha")
    stan_data_H0_alpha <- list(
      nM = length(M),
      M = M,
      X = X,
      alpha_gamma = alpha_gamma,
      beta_gamma  = beta_gamma,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      d_nu = d_nu,
      nu_lower = nu_lower
    )

    fit_H0_alpha <- rstan::stan(
      model_code = model_code_H0_alpha,
      data       = stan_data_H0_alpha,
      iter       = iter,
      warmup     = warmup,
      chains     = chains,
      seed       = seed,
      verbose    = verbose,
      refresh    = refresh,
      ...
    )

    print("Loading the Model: XM2Y_without_beta")
    stan_data_H0_beta <- list(
      nY = length(Y),
      Y = Y,
      M = M,
      X = X,
      alpha_gamma = alpha_gamma,
      beta_gamma  = beta_gamma,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      d_nu = d_nu,
      nu_lower = nu_lower
    )

    fit_H0_beta <- rstan::stan(
      model_code = model_code_H0_beta,
      data       = stan_data_H0_beta,
      iter       = iter,
      warmup     = warmup,
      chains     = chains,
      seed       = seed,
      verbose    = verbose,
      refresh    = refresh,
      ...
    )

    ## Posterior Samples of mediation effect
    samples_fit_H1_alpha <- extract(fit_H1_alpha)
    samples_fit_H1_beta <- extract(fit_H1_beta)
    samples_med <- samples_fit_H1_alpha$alpha * samples_fit_H1_beta$beta

    ## Bayes factor: mediation effect
    print("Bridge sampling started:")
    bridge_alpha_H0 <-  bridge_sampler(fit_H0_alpha,  repetitions = repetitions,
                                       method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                       maxiter = maxiter, silent = silent)
    bridge_alpha_H1 <-  bridge_sampler(fit_H1_alpha,  repetitions = repetitions,
                                       method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                       maxiter = maxiter, silent = silent)
    bf_alpha<- as.numeric(bf(bridge_alpha_H1,bridge_alpha_H0))[1]

    bridge_beta_H0 <-  bridge_sampler(fit_H0_beta,  repetitions = repetitions,
                                       method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                       maxiter = maxiter, silent = silent)
    bridge_beta_H1 <-  bridge_sampler(fit_H1_beta,  repetitions = repetitions,
                                       method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                       maxiter = maxiter, silent = silent)
    bf_beta<- as.numeric(bf(bridge_beta_H1,bridge_beta_H0))[1]

    bf_med <- (bf_alpha*bf_beta)/(q00+q10*bf_alpha+q01*bf_beta)

    ## Output
    output_list <- list(X2M_with_alpha = fit_H1_alpha,
                     XM2Y_with_beta = fit_H1_beta,
                     X2M_without_alpha = fit_H0_alpha,
                     XM2Y_without_beta = fit_H0_beta,
                     samples_med = samples_med,
                     bf_med = bf_med)
  }

  ################## Gamma-Only Model ########################
  if(model == "gamma-only"){
    model_code_H1_alpha <- "
functions {
  real log_likelihood_lpdf(vector m, vector x, real intercept, real alpha, real sigma, real gamma) {
    int n = num_elements(m);
    real c = sqrt(pi() / 2) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept - alpha * x[i] + sigma * c >= 0) {
        log_m_comp[i] = normal_lpdf((m[i] - intercept - alpha * x[i] + sigma * c) / (sigma * gamma) | 0, 1);
      } else {
        log_m_comp[i] = normal_lpdf((m[i] - intercept - alpha * x[i] + sigma * c) / (sigma * (1 / gamma)) | 0, 1);
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
  real alpha; // slope
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
}

model {
  // Prior: p(intercept, alpha, sigma)  1/sigma
  target += -log(sigma);

  // Prior: gamma
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
          - log_diff_exp(
              gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
              gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
            );

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, alpha, sigma, gamma);
}
    "

    model_code_H1_beta <- "
functions {
  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real beta, real tau, real sigma, real gamma) {
    int n = num_elements(y);
    real c = sqrt(pi() / 2) * (gamma - (1 / gamma));
    real log_l1 = n * log(2) - n * log(sigma) - n * log(gamma + (1 / gamma));
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] - beta * m[i] + sigma * c >= 0) {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i] - beta * m[i] + sigma * c) / (sigma * gamma) | 0, 1);
      } else {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i] - beta * m[i] + sigma * c) / (sigma * (1 / gamma)) | 0, 1);
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
  real beta; // slope
  real tau;
  real<lower=0> sigma; // standard deviation
  real<lower=gamma_lower, upper=gamma_upper> gamma;
}

model {
  // Prior: p(intercept, beta, tau, sigma)  1/sigma
  target += -log(sigma);

  // Prior: gamma
  target += gamma_lpdf(gamma | alpha_gamma, beta_gamma)
          - log_diff_exp(
              gamma_cdf(gamma_upper, alpha_gamma, beta_gamma),
              gamma_cdf(gamma_lower, alpha_gamma, beta_gamma)
            );

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, beta, tau, sigma, gamma);
}
    "

    model_code_H0_alpha <- "
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
  // Prior: p(intercept, sigma)  1/sigma
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
    "

    model_code_H0_beta <- "
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
  // Prior: p(intercept, tau, sigma)  1/sigma
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
    "

## H1: Mediation exists
print("Loading the Model: X2M_with_alpha")
stan_data_H1_alpha <- list(
  nM = length(M),
  M = M,
  X = X,
  alpha_gamma = alpha_gamma,
  beta_gamma  = beta_gamma,
  gamma_lower = gamma_lower,
  gamma_upper = gamma_upper
)

fit_H1_alpha <- rstan::stan(
  model_code = model_code_H1_alpha,
  data       = stan_data_H1_alpha,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

print("Loading the Model: XM2Y_with_beta")
stan_data_H1_beta <- list(
  nY = length(Y),
  Y = Y,
  M = M,
  X = X,
  alpha_gamma = alpha_gamma,
  beta_gamma  = beta_gamma,
  gamma_lower = gamma_lower,
  gamma_upper = gamma_upper
)

fit_H1_beta <- rstan::stan(
  model_code = model_code_H1_beta,
  data       = stan_data_H1_beta,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

## H0: Mediation does not exist
print("Loading the Model: X2M_without_alpha")
stan_data_H0_alpha <- list(
  nM = length(M),
  M = M,
  X = X,
  alpha_gamma = alpha_gamma,
  beta_gamma  = beta_gamma,
  gamma_lower = gamma_lower,
  gamma_upper = gamma_upper
)

fit_H0_alpha <- rstan::stan(
  model_code = model_code_H0_alpha,
  data       = stan_data_H0_alpha,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

print("Loading the Model: XM2Y_without_beta")
stan_data_H0_beta <- list(
  nY = length(Y),
  Y = Y,
  M = M,
  X = X,
  alpha_gamma = alpha_gamma,
  beta_gamma  = beta_gamma,
  gamma_lower = gamma_lower,
  gamma_upper = gamma_upper
)

fit_H0_beta <- rstan::stan(
  model_code = model_code_H0_beta,
  data       = stan_data_H0_beta,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

## Posterior Samples of mediation effect
samples_fit_H1_alpha <- extract(fit_H1_alpha)
samples_fit_H1_beta <- extract(fit_H1_beta)
samples_med <- samples_fit_H1_alpha$alpha * samples_fit_H1_beta$beta

## Bayes factor: mediation effect
print("Bridge sampling started:")
bridge_alpha_H0 <-  bridge_sampler(fit_H0_alpha,  repetitions = repetitions,
                                   method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                   maxiter = maxiter, silent = silent)
bridge_alpha_H1 <-  bridge_sampler(fit_H1_alpha,  repetitions = repetitions,
                                   method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                   maxiter = maxiter, silent = silent)
bf_alpha<- as.numeric(bf(bridge_alpha_H1,bridge_alpha_H0))[1]

bridge_beta_H0 <-  bridge_sampler(fit_H0_beta,  repetitions = repetitions,
                                  method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                  maxiter = maxiter, silent = silent)
bridge_beta_H1 <-  bridge_sampler(fit_H1_beta,  repetitions = repetitions,
                                  method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                  maxiter = maxiter, silent = silent)
bf_beta<- as.numeric(bf(bridge_beta_H1,bridge_beta_H0))[1]

bf_med <- (bf_alpha*bf_beta)/(q00+q10*bf_alpha+q01*bf_beta)

## Output
output_list <- list(X2M_with_alpha = fit_H1_alpha,
                    XM2Y_with_beta = fit_H1_beta,
                    X2M_without_alpha = fit_H0_alpha,
                    XM2Y_without_beta = fit_H0_beta,
                    samples_med = samples_med,
                    bf_med = bf_med)

  }

  ################## nu-Only Model ########################
  if(model == "nu-only"){
    model_code_H1_alpha <- "
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
  // Prior: p(intercept, alpha, sigma)  1/sigma
  target += -log(sigma);

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, alpha, sigma, nu);
}
    "

    model_code_H1_beta <- "
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
  // Prior: p(intercept, beta, tau, sigma)  1/sigma
  target += -log(sigma);

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, beta, tau, sigma, nu);
}
    "

    model_code_H0_alpha <- "
functions {
  real log_likelihood_lpdf(vector m, vector x, real intercept, real sigma, real nu) {
    int n = num_elements(m);
    real log_l1 = - n * log(sigma);
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept >= 0) {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept) / (sigma) | nu, 0, 1);
      } else {
        log_m_comp[i] = student_t_lpdf((m[i] - intercept) / (sigma) | nu, 0, 1);
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
  real<lower=0> sigma; // standard deviation
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, sigma)  1/sigma
  target += -log(sigma);

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, sigma, nu);
}
    "

    model_code_H0_beta <- "
functions {
  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real tau, real sigma, real nu) {
    int n = num_elements(y);
    real log_l1 = - n * log(sigma);
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] >= 0) {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i]) / (sigma) | nu, 0, 1);
      } else {
        log_y_comp[i] = student_t_lpdf((y[i] - intercept - tau * x[i]) / (sigma) | nu, 0, 1);
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
  real tau;
  real<lower=0> sigma; // standard deviation
  real<lower=nu_lower> nu;
}

model {
  // Prior: p(intercept, tau, sigma)  1/sigma
  target += -log(sigma);

  // Prior: nu
  target += exponential_lpdf(nu | d_nu)
          - exponential_lccdf(nu_lower | d_nu);

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, tau, sigma, nu);
}
    "

## H1: Mediation exists
print("Loading the Model: X2M_with_alpha")
stan_data_H1_alpha <- list(
  nM = length(M),
  M = M,
  X = X,
  d_nu = d_nu,
  nu_lower = nu_lower
)

fit_H1_alpha <- rstan::stan(
  model_code = model_code_H1_alpha,
  data       = stan_data_H1_alpha,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

print("Loading the Model: XM2Y_with_beta")
stan_data_H1_beta <- list(
  nY = length(Y),
  Y = Y,
  M = M,
  X = X,
  d_nu = d_nu,
  nu_lower = nu_lower
)

fit_H1_beta <- rstan::stan(
  model_code = model_code_H1_beta,
  data       = stan_data_H1_beta,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

## H0: Mediation does not exist
print("Loading the Model: X2M_without_alpha")
stan_data_H0_alpha <- list(
  nM = length(M),
  M = M,
  X = X,
  d_nu = d_nu,
  nu_lower = nu_lower
)

fit_H0_alpha <- rstan::stan(
  model_code = model_code_H0_alpha,
  data       = stan_data_H0_alpha,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

print("Loading the Model: XM2Y_without_beta")
stan_data_H0_beta <- list(
  nY = length(Y),
  Y = Y,
  M = M,
  X = X,
  d_nu = d_nu,
  nu_lower = nu_lower
)

fit_H0_beta <- rstan::stan(
  model_code = model_code_H0_beta,
  data       = stan_data_H0_beta,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

## Posterior Samples of mediation effect
samples_fit_H1_alpha <- extract(fit_H1_alpha)
samples_fit_H1_beta <- extract(fit_H1_beta)
samples_med <- samples_fit_H1_alpha$alpha * samples_fit_H1_beta$beta

## Bayes factor: mediation effect
print("Bridge sampling started:")
bridge_alpha_H0 <-  bridge_sampler(fit_H0_alpha,  repetitions = repetitions,
                                   method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                   maxiter = maxiter, silent = silent)
bridge_alpha_H1 <-  bridge_sampler(fit_H1_alpha,  repetitions = repetitions,
                                   method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                   maxiter = maxiter, silent = silent)
bf_alpha<- as.numeric(bf(bridge_alpha_H1,bridge_alpha_H0))[1]

bridge_beta_H0 <-  bridge_sampler(fit_H0_beta,  repetitions = repetitions,
                                  method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                  maxiter = maxiter, silent = silent)
bridge_beta_H1 <-  bridge_sampler(fit_H1_beta,  repetitions = repetitions,
                                  method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                  maxiter = maxiter, silent = silent)
bf_beta<- as.numeric(bf(bridge_beta_H1,bridge_beta_H0))[1]

bf_med <- (bf_alpha*bf_beta)/(q00+q10*bf_alpha+q01*bf_beta)

## Output
output_list <- list(X2M_with_alpha = fit_H1_alpha,
                    XM2Y_with_beta = fit_H1_beta,
                    X2M_without_alpha = fit_H0_alpha,
                    XM2Y_without_beta = fit_H0_beta,
                    samples_med = samples_med,
                    bf_med = bf_med)
  }


  ################## Normal Model ########################
  if(model == "normal"){
    model_code_H1_alpha <- "
functions {
  real log_likelihood_lpdf(vector m, vector x, real intercept, real alpha, real sigma) {
    int n = num_elements(m);
    real log_l1 = - n * log(sigma);
    vector[n] log_m_comp;

    for (i in 1:n) {
      if (m[i] - intercept - alpha * x[i]>= 0) {
        log_m_comp[i] = normal_lpdf((m[i] - intercept - alpha * x[i]) / (sigma) | 0, 1);
      } else {
        log_m_comp[i] = normal_lpdf((m[i] - intercept - alpha * x[i]) / (sigma) | 0, 1);
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
  real alpha; // slope
  real<lower=0> sigma; // standard deviation
}

model {
  // Prior: p(intercept, alpha, sigma)  1/sigma
  target += -log(sigma);

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, alpha, sigma);
}
    "

    model_code_H1_beta <- "
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
  // Prior: p(intercept, beta, tau, sigma)  1/sigma
  target += -log(sigma);

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, beta, tau, sigma);
}
    "

    model_code_H0_alpha <- "
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
  // Prior: p(intercept, sigma)  1/sigma
  target += -log(sigma);

  // Likelihood
  target += log_likelihood_lpdf(M | X, intercept, sigma);
}
    "

    model_code_H0_beta <- "
functions {
  real log_likelihood_lpdf(vector y, vector m, vector x, real intercept, real tau, real sigma) {
    int n = num_elements(y);
    real log_l1 = - n * log(sigma);
    vector[n] log_y_comp;

    for (i in 1:n) {
      if (y[i] - intercept - tau * x[i] >= 0) {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i]) / (sigma) | 0, 1);
      } else {
        log_y_comp[i] = normal_lpdf((y[i] - intercept - tau * x[i]) / (sigma) | 0, 1);
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
  real tau;
  real<lower=0> sigma; // standard deviation
}

model {
  // Prior: p(intercept, tau, sigma)  1/sigma
  target += -log(sigma);

  // Likelihood
  target += log_likelihood_lpdf(Y | M, X, intercept, tau, sigma);
}
    "

## H1: Mediation exists
print("Loading the Model: X2M_with_alpha")
stan_data_H1_alpha <- list(
  nM = length(M),
  M = M,
  X = X
)

fit_H1_alpha <- rstan::stan(
  model_code = model_code_H1_alpha,
  data       = stan_data_H1_alpha,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

print("Loading the Model: XM2Y_with_beta")
stan_data_H1_beta <- list(
  nY = length(Y),
  Y = Y,
  M = M,
  X = X
)

fit_H1_beta <- rstan::stan(
  model_code = model_code_H1_beta,
  data       = stan_data_H1_beta,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

## H0: Mediation does not exist
print("Loading the Model: X2M_without_alpha")
stan_data_H0_alpha <- list(
  nM = length(M),
  M = M,
  X = X
)

fit_H0_alpha <- rstan::stan(
  model_code = model_code_H0_alpha,
  data       = stan_data_H0_alpha,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

print("Loading the Model: XM2Y_without_beta")
stan_data_H0_beta <- list(
  nY = length(Y),
  Y = Y,
  M = M,
  X = X
)

fit_H0_beta <- rstan::stan(
  model_code = model_code_H0_beta,
  data       = stan_data_H0_beta,
  iter       = iter,
  warmup     = warmup,
  chains     = chains,
  seed       = seed,
  verbose    = verbose,
  refresh    = refresh,
  ...
)

## Posterior Samples of mediation effect
samples_fit_H1_alpha <- extract(fit_H1_alpha)
samples_fit_H1_beta <- extract(fit_H1_beta)
samples_med <- samples_fit_H1_alpha$alpha * samples_fit_H1_beta$beta

## Bayes factor: mediation effect
print("Bridge sampling started:")
bridge_alpha_H0 <-  bridge_sampler(fit_H0_alpha,  repetitions = repetitions,
                                   method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                   maxiter = maxiter, silent = silent)
bridge_alpha_H1 <-  bridge_sampler(fit_H1_alpha,  repetitions = repetitions,
                                   method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                   maxiter = maxiter, silent = silent)
bf_alpha<- as.numeric(bf(bridge_alpha_H1,bridge_alpha_H0))[1]

bridge_beta_H0 <-  bridge_sampler(fit_H0_beta,  repetitions = repetitions,
                                  method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                  maxiter = maxiter, silent = silent)
bridge_beta_H1 <-  bridge_sampler(fit_H1_beta,  repetitions = repetitions,
                                  method = method, cores = bridgesampling_cores, use_neff = use_neff,
                                  maxiter = maxiter, silent = silent)
bf_beta<- as.numeric(bf(bridge_beta_H1,bridge_beta_H0))[1]

bf_med <- (bf_alpha*bf_beta)/(q00+q10*bf_alpha+q01*bf_beta)

## Output
output_list <- list(X2M_with_alpha = fit_H1_alpha,
                    XM2Y_with_beta = fit_H1_beta,
                    X2M_without_alpha = fit_H0_alpha,
                    XM2Y_without_beta = fit_H0_beta,
                    samples_med = samples_med,
                    bf_med = bf_med)

  }

  return(output_list)

}
