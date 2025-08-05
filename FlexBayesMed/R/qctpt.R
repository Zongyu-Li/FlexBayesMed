#' Quantile Function of the Centred Two-Piece Student \emph{t}-Distribution (CTPT)
#'
#' Inverts the CDF of the CTPT distribution to find \eqn{\varepsilon} such that
#' \eqn{P(X \le \varepsilon) = p}.
#'
#' The distribution is defined by a piecewise CDF:
#' \deqn{
#'   F(\varepsilon|\gamma,\nu) =
#'   \begin{cases}
#'     \displaystyle
#'     \frac{2}{1 + \gamma^2}\,F_\nu\bigl[\gamma\,(\varepsilon + m)\bigr], & \varepsilon < -m,\\[2em]
#'     \displaystyle
#'     \frac{2\,\gamma^2}{1 + \gamma^2}\,F_\nu\!\Bigl(\tfrac{\varepsilon + m}{\gamma}\Bigr)
#'       \;+\;\frac{1 - \gamma^2}{1 + \gamma^2}, & \varepsilon \ge -m,
#'   \end{cases}
#' }
#' where \eqn{m = m(\gamma,\nu)} is the mean-adjustment term, and \eqn{F_\nu} is the
#' standard Student t (if \eqn{\nu > 1}) or standard normal (if \code{nu="normal"}) CDF.
#'
#' @param p Numeric vector of probabilities in \eqn{[0,1]}.
#' @param gamma Positive scalar (\eqn{\gamma > 0}). The skewness parameter.
#' @param nu Either a numeric \eqn{>1} (degrees of freedom) or the string \code{"normal"}.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \le q)}.
#' @param log.p Logical; if \code{TRUE}, \code{p} is treated as \eqn{\log(p)}.
#'
#' @return Numeric vector of quantiles corresponding to the given probabilities.
#'
#' @seealso \code{\link{pctpt}} for the CDF and \code{\link{dctpt}} for the PDF.
#'
#' @examples
#' qctpt(p = 0.975, gamma = 1, nu = "normal") # Returns 1.96
#' qctpt(p = 0.025, gamma = 2, nu = 2) # Returns -3.40
#'
#' @export
#' @importFrom stats qnorm qt
qctpt <- function(p, gamma, nu, lower.tail = TRUE, log.p = FALSE) {
  if (!is.numeric(gamma) || gamma <= 0) {
    stop("gamma must be a positive numeric value.")
  }
  if (is.character(nu)) {
    if (nu != "normal") {
      stop("nu must be 'normal' or a numeric value > 1.")
    }
  } else {
    if (nu <= 1) {
      stop("For Student t, 'nu' must be > 1.")
    }
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("All probabilities 'p' must be in [0, 1].")
  }

  m <- mean_two_piece_student(gamma, nu)

  # Boundary probability:  p_boundary = F_X(-m) = 1 / (1 + gamma^2)
  p_boundary <- 1 / (1 + gamma^2)

  # Inverse distribution function for standard t or normal
  base_inv <- if (is.character(nu) && nu == "normal") {
    function(prob) qnorm(prob)
  } else {
    function(prob) qt(prob, df = nu)
  }

  q_out <- numeric(length(p))

  # Left piece
  idx_left <- (p <= p_boundary)
  p_left   <- p[idx_left]
  if (any(idx_left)) {
    arg_left <- (p_left * (1 + gamma^2)) / 2
    q_out[idx_left] <- (1 / gamma) * base_inv(arg_left) - m
  }

  # Right piece
  idx_right <- !idx_left
  p_right   <- p[idx_right]
  if (any(idx_right)) {
    numerator <- p_right * (1 + gamma^2) - (1 - gamma^2)
    denominator <- 2 * gamma^2
    arg_right <- numerator / denominator
    q_out[idx_right] <- gamma * base_inv(arg_right) - m
  }

  return(q_out)
}
