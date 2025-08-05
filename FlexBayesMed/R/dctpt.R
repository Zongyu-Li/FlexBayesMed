#' Density of the Centered Two-Piece Student \eqn{t}-Distribution (CTPT)
#'
#' @param x Numeric vector. Quantiles at which to compute the density.
#' @param gamma Positive numeric scalar. Skewness parameter.
#' @param nu Numeric scalar \eqn{\nu > 1} (degrees of freedom) or "normal".
#' @param log Logical. Return log-density? Default: `FALSE`.
#'
#' @return Numeric vector of densities.
#'
#' @details
#' The CTPT distribution is defined by shifting the original two-piece Student \eqn{t}-distribution
#' (with skewness parameter \eqn{\gamma} and degrees of freedom \eqn{\nu}) by its mean, ensuring the
#' resulting distribution has zero mean. The density of the CTPT distribution is given by:
#' \deqn{%
#'   \begin{aligned}
#'       p(\varepsilon|\gamma,\nu) = \frac{2}{\gamma + \frac{1}{\gamma}} \Bigg[& f_\nu\left(\frac{\varepsilon + m(\gamma,\nu)}{\gamma}\right) I_{[-m(\gamma,\nu),+\infty)}(\varepsilon) \;+ \\
#'       &f_\nu\left(\gamma\big(\varepsilon + m(\gamma,\nu)\big)\right) I_{(-\infty,-m(\gamma,\nu))}(\varepsilon) \Bigg], \quad \gamma \in \mathbb{R}_+, \; \nu > 1,
#'   \end{aligned}
#' }
#' where:
#' \itemize{
#'   \item \eqn{f_\nu} is the Student \eqn{t}-density with \eqn{\nu} degrees of freedom.
#'   \item \eqn{m(\gamma,\nu)} is the mean of the unadjusted two-piece distribution, computed by \code{\link{mean_two_piece_student}}.
#'   \item \eqn{I_A(\varepsilon)} denotes the indicator function for set \eqn{A}.
#' }
#'
#' @examples
#' x <- seq(-5, 5, length.out = 100)
#' density <- dctpt(x, gamma = 2, nu = 5)
#' plot(x, density, type = "l")
#' @importFrom stats dnorm dt
#' @export
dctpt <- function(x, gamma, nu, log = FALSE) {
  # Parameter checks
  if (!is.numeric(gamma) || gamma <= 0)
    stop("gamma must be a positive numeric value.")
  if (is.character(nu)) {
    if (nu != "normal")
      stop("nu must be 'normal' or a numeric value > 1.")
  } else if (nu <= 1) {
    stop("nu must be > 1 or 'normal'.")
  }

  m <- mean_two_piece_student(gamma, nu)
  y <- x + m  # Shift to original two-piece variable

  # Compute density based on distribution type
  if (is.character(nu)) {
    # Normal case
    left_part <- y < 0
    density <- numeric(length(y))
    density[left_part] <- (2 / (gamma + 1/gamma)) * dnorm(gamma * y[left_part])
    density[!left_part] <- (2 / (gamma + 1/gamma)) * dnorm(y[!left_part] / gamma)
  } else {
    # Student t case
    left_part <- y < 0
    density <- numeric(length(y))
    density[left_part] <- (2 / (gamma + 1/gamma)) * dt(gamma * y[left_part], df = nu)
    density[!left_part] <- (2 / (gamma + 1/gamma)) * dt(y[!left_part] / gamma, df = nu)
  }

  if (log) density <- log(density)
  return(density)
}
