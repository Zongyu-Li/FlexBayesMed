#' Variance of the Uncentered Two-Piece Student's t (or Normal) Distribution
#'
#' This function calculates the variance of the unadjusted two-piece Student's t distribution
#' (see Proposition 2 in ...) for \eqn{\gamma > 0} and \eqn{\nu > 2}.
#' It also handles the limiting "normal" case (i.e., \eqn{\nu \to \infty}).
#'
#' @param gamma A positive numeric value giving the skewness parameter (must be \eqn{\gamma > 0}).
#' @param nu Either:
#'   \itemize{
#'     \item A numeric value \code{> 2}, giving the degrees of freedom \eqn{\nu} for the Student's t
#'       (required for a finite variance).
#'     \item The character string \code{"normal"} to indicate the limiting two-piece normal case
#'       (i.e., \eqn{\nu = \infty}).
#'   }
#'
#' @details
#' When \code{nu} is a finite numeric value \eqn{> 2}, the variance formula is:
#' \deqn{
#'   \mathrm{Var}(\varepsilon \mid \gamma,\nu)
#'   =
#'     \frac{\nu}{\nu - 2}
#'     \Bigl(\gamma^2 - 1 + \tfrac{1}{\gamma^2}\Bigr)
#'     \;-\;
#'     m(\gamma,\nu)^2,
#' }
#' where
#' \deqn{
#'   m(\gamma,\nu)
#'   =
#'     \frac{2\,\nu\,\Gamma\bigl((\nu + 1)/2\bigr)}{\sqrt{\pi\,\nu}\,(\nu-1)\,\Gamma\bigl(\nu/2\bigr)}
#'     \,\bigl(\gamma - \tfrac{1}{\gamma}\bigr).
#' }
#'
#' When \code{nu = "normal"}, the variance is that of the unadjusted two-piece normal distribution:
#' \deqn{
#'   \mathrm{Var}(\varepsilon \mid \gamma,\nu=\infty)
#'   =
#'     \bigl(\gamma^2 - 1 + \tfrac{1}{\gamma^2}\bigr)
#'     \;-\;
#'     \Bigl(\sqrt{\tfrac{2}{\pi}}\bigl(\gamma - \tfrac{1}{\gamma}\bigr)\Bigr)^2.
#' }
#'
#' @return A numeric value giving the variance of the unadjusted two-piece distribution.
#'
#' @examples
#' # Variance for gamma=2, nu=5
#' var_two_piece_student(2, 5)
#'
#' # Variance in the limiting normal case for gamma=2
#' var_two_piece_student(2, "normal")
#'
#' @export

var_two_piece_student <- function(gamma, nu) {

  if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("gamma must be a single positive numeric value.")
  }

  if (identical(nu, "normal")) {
    m_norm <- sqrt(2 / pi) * (gamma - 1/gamma)
    A2_norm <- gamma^2 - 1 + 1/gamma^2
    var_val <- A2_norm - m_norm^2

  } else {
    if (!is.numeric(nu) || length(nu) != 1) {
      stop("nu must be either 'normal' or a single numeric value > 2.")
    }
    nu_val <- as.numeric(nu)
    if (nu_val <= 2) {
      stop("nu must be strictly greater than 2 for variance to exist.")
    }

    m_num <- 2 * nu_val * gamma((nu_val + 1) / 2)
    m_den <- sqrt(pi * nu_val) * (nu_val - 1) * gamma(nu_val / 2)
    m_val <- (m_num / m_den) * (gamma - 1/gamma)

    # Use Proposition 2 (Li et al., 2025)
    var_val <- (nu_val / (nu_val - 2)) * (gamma^2 - 1 + 1/gamma^2) - m_val^2
  }

  return(var_val)
}
