#' Generate Samples from the Mean-Adjusted Two-Piece Student's \eqn{t}-Distribution (MATPT)
#'
#' Produces independent and identically distributed (i.i.d.) samples from the MATPT distribution,
#' which is a zero-mean adjusted version of the two-piece Student's \eqn{t}-distribution (or normal).
#' The adjustment ensures that the error terms satisfy the zero-mean assumption, which is needed in cases such as mean regression or mediation analysis.
#'
#' @param N Positive integer. The number of samples to generate.
#' @param gamma Positive numeric scalar. Skewness parameter (\eqn{\gamma > 0}).
#' @param nu Either a numeric scalar \eqn{\nu > 1} (degrees of freedom) or the string "normal" for
#'           the normal limit (\eqn{\nu \to \infty}).
#'
#' @return A numeric vector of length \eqn{N} containing the samples.
#'
#' @details
#' The MATPT distribution is defined by shifting the original two-piece Student's \eqn{t}-distribution
#' (with skewness parameter \eqn{\gamma} and degrees of freedom \eqn{\nu}) by its mean, ensuring the
#' resulting distribution has zero mean. The density of the MATPT distribution is given by:
#' \deqn{%
#'   \begin{aligned}
#'       p(\varepsilon|\gamma,\nu) = \frac{2}{\gamma + \frac{1}{\gamma}} \Bigg[& f_\nu\left(\frac{\varepsilon + m(\gamma,\nu)}{\gamma}\right) I_{[-m(\gamma,\nu),+\infty)}(\varepsilon) \;+ \\
#'       &f_\nu\left(\gamma\big(\varepsilon + m(\gamma,\nu)\big)\right) I_{(-\infty,-m(\gamma,\nu))}(\varepsilon) \Bigg], \quad \gamma \in \mathbb{R}_+, \; \nu > 1,
#'   \end{aligned}
#' }
#' where:
#' \itemize{
#'   \item \eqn{f_\nu} is the Student's \eqn{t}-density with \eqn{\nu} degrees of freedom.
#'   \item \eqn{m(\gamma,\nu)} is the mean of the unadjusted two-piece distribution, computed by \code{\link{mean_two_piece_student}}.
#'   \item \eqn{I_A(\varepsilon)} denotes the indicator function for set \eqn{A}.
#' }
#'
#' @references
#' TBA
#'
#' @examples
#' # Generate samples from MATPT with gamma=2, nu=5 (Student's t case)
#' samples_t <- rmatpt(N = 1000, gamma = 2, nu = 5)
#' hist(samples_t, breaks = 30, main = "MATPT (nu=5)")
#'
#' # Generate samples from MATPT with gamma=2 (normal case)
#' samples_norm <- rmatpt(N = 1000, gamma = 2, nu = "normal")
#' hist(samples_norm, breaks = 30, main = "MATPT (normal)")
#'
#' @seealso
#' \code{\link{mean_two_piece_student}} for computing the mean of the unadjusted distribution.
#'
#' @importFrom stats runif rnorm rt
#' @export

rmatpt <- function(N, gamma, nu) {
  if (!is.numeric(gamma) || gamma <= 0) {
    stop("'gamma' must be a positive numeric value.")
  }

  ## If nu is character == "normal", use the normal sampling
  if (is.character(nu)) {
    if (nu == "normal") {
      m <- mean_two_piece_student(gamma,nu="normal")
      p_left <- 1/(gamma^2 + 1)
      x <- numeric(N)

      for (i in seq_len(N)) {
        u <- runif(1)
        if (u <= p_left) {
          repeat {
            a <- rnorm(1)
            if (a < 0) {
              a <- a/gamma
              break
            }
          }
        } else {
          repeat {
            a <- rnorm(1)
            if (a >= 0) {
              a <- a*gamma
              break
            }
          }
        }
        x[i] <- a
      }
      x <- x - m
      return(x)

    } else {
      stop("'nu' must be a positive number or the string 'normal'.")
    }
    # when nu < \infty, use the Student t sampling
  } else if (is.numeric(nu)) {
    if (nu <= 0) {
      stop("'nu' must be a positive numeric value (or 'normal').")
    }
    m <- mean_two_piece_student(gamma, nu)
    p_left <- 1/(gamma^2 + 1)
    x <- numeric(N)

    for (i in seq_len(N)) {
      u <- runif(1)
      if (u <= p_left) {
        repeat {
          a <- rt(1, df=nu)
          if (a < 0) {
            a <- a/gamma
            break
          }
        }
      } else {
        repeat {
          a <- rt(1, df=nu)
          if (a >= 0) {
            a <- a*gamma
            break
          }
        }
      }
      x[i] <- a
    }
    x <- x - m
    return(x)
  } else {
    stop("'nu' must be numeric or the string 'normal'.")
  }
}
