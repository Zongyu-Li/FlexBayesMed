#' Cumulative Distribution Function (CDF) of the Centred Two-Piece Student \emph{t}-Distribution (CTPT)
#'
#'
#' @param q Numeric vector of quantiles.
#' @param gamma Positive scalar. Skewness parameter (\eqn{\gamma > 0}).
#' @param nu Degrees of freedom. Either numeric \code{>1} or \code{"normal"} for the Gaussian limit.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \le q)}.
#' @param log.p Logical; if \code{TRUE}, returns log probabilities.
#'
#' @section Details:
#' The CDF is piecewise-defined with split point at \eqn{\varepsilon = -m(\gamma,\nu)}:
#'
#' For \eqn{\varepsilon < -m(\gamma,\nu)}:
#' \deqn{
#'   F(\varepsilon \mid \gamma,\nu)
#'   = \frac{2}{1+\gamma^2}\, F_\nu\!\bigl[\gamma\,(\varepsilon + m(\gamma,\nu))\bigr]
#' }
#'
#' For \eqn{\varepsilon \ge -m(\gamma,\nu)}:
#' \deqn{
#'   F(\varepsilon \mid \gamma,\nu)
#'   = \frac{2\gamma^2}{1+\gamma^2}\, F_\nu\!\bigl[\tfrac{\varepsilon + m(\gamma,\nu)}{\gamma}\bigr]
#'     \;+\; \frac{1-\gamma^2}{1+\gamma^2}
#' }
#'
#' Where:
#' \itemize{
#'   \item \eqn{m(\gamma,\nu)} is the mean-adjustment term (see Proposition 1 in manuscript)
#'   \item \eqn{F_\nu} = CDF of Student t (\eqn{\nu} df) or standard normal (\code{nu = "normal"})
#'   \item Gaussian mean: \eqn{m(\gamma) = \sqrt{\tfrac{2}{\pi}}\,\bigl(\gamma - \tfrac{1}{\gamma}\bigr)}
#' }
#'
#' @references TBA
#'
#' @seealso
#' [mean_two_piece_student()] for the mean-adjustment computation
#'
#' @examples
#' # Symmetric case (gamma=1) matches standard t
#' pctpt(q = 0, gamma = 1, nu = 5)  # Returns 0.5
#' pt(0, df = 5)                      # Returns 0.5
#'
#' # Skewed case
#' pctpt(q = 0, gamma = 2, nu = 5)
#'
#' # Normal limit
#' pctpt(q = 0, gamma = 1.5, nu = "normal")
#'
#' @export
#' @importFrom stats pnorm pt


pctpt <- function(q, gamma, nu, lower.tail = TRUE, log.p = FALSE) {
  # Parameter validation
  if (!is.numeric(gamma) || gamma <= 0)
    stop("gamma must be a positive numeric value.")

  if (is.character(nu)) {
    if (nu != "normal")
      stop("nu must be 'normal' or a numeric value > 1.")
  } else if (nu <= 1) {
    stop("For Student t sampling, nu must be > 1.")
  }

  # Compute mean adjustment
  m <- mean_two_piece_student(gamma, nu)

  # Shift quantiles
  y <- q + m

  # Select base distribution
  base_cdf <- if (is.character(nu)) {
    function(z) pnorm(z, lower.tail = TRUE)
  } else {
    function(z) pt(z, df = nu, lower.tail = TRUE)
  }

  # Compute CDF
  cdf <- numeric(length(q))
  left <- (y < 0)
  right <- !left

  cdf[left] <- (2 / (1 + gamma^2)) * base_cdf(gamma * y[left])
  cdf[right] <- (2 * gamma^2 / (1 + gamma^2)) * base_cdf(y[right] / gamma) +
    ((1 - gamma^2) / (1 + gamma^2))

  # Adjust for tail and log
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)

  return(cdf)
}
