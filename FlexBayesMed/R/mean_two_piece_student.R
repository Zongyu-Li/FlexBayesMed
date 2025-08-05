#' Mean of the Uncentered Two-Piece Student \emph{t} (or Normal) Distribution
#'
#' Computes the theoretical mean
#' \eqn{m(\gamma,\nu)} of the two-piece Student
#' \eqn{t_\nu} distribution with skewness parameter
#' \eqn{\gamma}. If `nu = "normal"` the function returns the
#' two-piece \emph{normal} limit
#' \eqn{m^{\textrm{normal}}(\gamma)=\sqrt{2/\pi}\,(\gamma-1/\gamma)} as \eqn{\nu \to \infty}.
#' This corresponds to Eq. (1) in Fernández & Steel (1998).
#'
#' @name mean_two_piece_student
#' @param gamma Positive numeric scalar. Skewness parameter \eqn{\gamma>0}.
#' @param nu    Either a numeric scalar \eqn{\nu>1} (degrees of freedom)
#'              \emph{or} the literal string `"normal"` for the
#'              \eqn{\nu \to \infty} limit.
#'
#' @return A numeric scalar giving the mean \eqn{m(\gamma,\nu)}.
#'
#' @details The closed-form expression for finite \eqn{\nu} is
#' \deqn{
#'   m(\gamma,\nu)=
#'     \frac{2\nu\Gamma\!\left(\frac{\nu+1}{2}\right)}
#'          {\sqrt{\pi\nu}\,(\nu-1)\Gamma\!\left(\frac{\nu}{2}\right)}
#'     \left(\gamma-\frac{1}{\gamma}\right).
#' }
#' For the normal limit (`nu = "normal"`) the mean simplifies to
#' \eqn{\sqrt{2/\pi}\,(\gamma-1/\gamma)}.
#'
#' @references
#' Fernández, C. & Steel, M. F. J. (1998).
#' “On Bayesian Modeling of Fat Tails and Skewness.”
#' \emph{Journal of the American Statistical Association},
#' \strong{93}(441), 359–371.
#'
#' @examples
#' mean_two_piece_student(gamma = 2, nu = 5)        # finite-nu case, approximately 1.424
#' mean_two_piece_student(gamma = 3, nu = "normal") # normal limit, approximately 2.128
#'
#' @importFrom Rmpfr mpfr
#' @export

mean_two_piece_student <- function(gamma, nu){

  if (!is.numeric(gamma) || length(gamma) != 1L || gamma <= 0)
    stop("'gamma' must be a positive numeric scalar.", call. = FALSE)

  if (is.character(nu) && identical(nu, "normal")) {
    return(sqrt(2 / pi) * (gamma - 1 / gamma))
  }

  if (!is.numeric(nu) || length(nu) != 1L || nu <= 1)
    stop("'nu' must be > 1 (or the string \"normal\").", call. = FALSE)

  nu_mpfr   <- Rmpfr::mpfr(nu, precBits = 512)
  log_num   <- log(2) + log(nu_mpfr) + lgamma((nu_mpfr + 1) / 2)
  log_den   <- 0.5 * (log(pi) + log(nu_mpfr)) + log(nu_mpfr - 1) + lgamma(nu_mpfr / 2)
  log_ratio <- log_num - log_den

  return(as.numeric(exp(log_ratio) * (gamma - 1 / gamma)))
}
