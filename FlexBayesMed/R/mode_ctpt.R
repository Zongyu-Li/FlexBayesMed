#' Mode of the Centered Two-Piece Student \emph{t} (or Normal) Distribution
#'
#' Computes the theoretical mode
#' \eqn{mode(\gamma,\nu)} of the centred two-piece Student
#' \eqn{t_\nu} distribution with skewness parameter
#' \eqn{\gamma}. If `nu = "normal"` the function returns the
#' two-piece \emph{normal} limit
#' \eqn{mode^{\textrm{normal}}(\gamma)=-\sqrt{2/\pi}\,(\gamma-1/\gamma)} as \eqn{\nu \to \infty}.
#'
#' @name mode_ctpt
#' @param gamma Positive numeric scalar. Skewness parameter \eqn{\gamma>0}.
#' @param nu    Either a numeric scalar \eqn{\nu>1} (degrees of freedom)
#'              \emph{or} the literal string `"normal"` for the
#'              \eqn{\nu \to \infty} limit.
#'
#' @return A numeric scalar giving the mean \eqn{mode(\gamma,\nu)}.
#'
#' @details The closed-form expression for finite \eqn{\nu} is
#' \deqn{
#'   mode(\gamma,\nu)= -
#'     \frac{2\nu\Gamma\!\left(\frac{\nu+1}{2}\right)}
#'          {\sqrt{\pi\nu}\,(\nu-1)\Gamma\!\left(\frac{\nu}{2}\right)}
#'     \left(\gamma-\frac{1}{\gamma}\right).
#' }
#' For the normal limit (`nu = "normal"`) the mean simplifies to
#' \eqn{-\sqrt{2/\pi}\,(\gamma-1/\gamma)}.
#'
#'
#' @examples
#' mode_ctpt(gamma = 2, nu = 5)        # finite-nu case, approximately -1.424
#' mode_ctpt(gamma = 3, nu = "normal") # normal limit, approximately -2.128
#'
#' @importFrom Rmpfr mpfr
#' @export

mode_ctpt <- function(gamma, nu){
  return(-mean_two_piece_student(gamma,nu))
}
