#' Calculate the Fisher's moment coefficient of skewness for CTPT
#'
#' This function calculates the skewness of the Mean-Adjusted Two-Piece Student's t (CTPT) distribution
#' for a given skewness parameter \code{gamma > 0} and degrees of freedom \code{nu > 3}.
#' It also handles the limiting normal case if \code{nu = "normal"}.
#'
#' @param gamma A positive numeric value giving the skewness parameter.
#' @param nu Either:
#'   \itemize{
#'     \item A numeric value \code{> 3}, giving the degrees of freedom for the Student's t
#'       (required for a finite skewness).
#'     \item The character string \code{"normal"} to indicate the limiting two-piece normal case.
#'   }
#'
#' @details
#' The formula for finite \code{nu} follows Proposition 3 in your manuscript.
#' When \code{nu = "normal"}, the skewness is computed under the limiting two-piece normal case.
#'
#' @return A numeric value giving the Fisher's moment coefficient of skewness.
#'
#' @examples
#' # Skewness for gamma=2, nu=5
#' skewness_ctpt(gamma = 2, nu = 5)
#'
#' # Skewness in the normal limit for gamma=2
#' skewness_ctpt(gamma = 2, nu = "normal")
#'
#' @export
skewness_ctpt <- function(gamma, nu = "normal") {

  # Check gamma
  if (gamma <= 0) stop("gamma must be positive.")

  # Check if nu is the "normal" string
  if (identical(nu, "normal")) {
    # --- Normal-limit case ---

    # Mean of unadjusted two-piece normal
    m_norm <- sqrt(2 / pi) * (gamma - 1/gamma)
    # 2nd raw moment
    A2_norm <- gamma^2 - 1 + 1/gamma^2
    # 3rd raw moment
    A3_norm <- 2 * sqrt(2 / pi) * ((gamma^4 - gamma^(-4)) / (gamma + 1/gamma))

    # Convert raw moments to central moments
    EX2 <- A2_norm - m_norm^2
    EX3 <- A3_norm - 3 * m_norm * A2_norm + 2 * m_norm^3

    return(EX3 / (EX2^(1.5)))

  } else {
    # Otherwise, we expect nu to be numeric
    if (!is.numeric(nu)) {
      stop("nu must be either 'normal' or a numeric value > 3.")
    }

    # Convert to numeric if not already
    nu_val <- as.numeric(nu)

    # Check degrees of freedom
    if (nu_val <= 3) {
      stop("For finite skewness, nu must be strictly greater than 3.")
    }

    # --- Finite nu case ---
    # Mean of unadjusted two-piece t
    m_num <- 2 * nu_val * gamma((nu_val + 1) / 2)
    m_den <- sqrt(pi * nu_val) * (nu_val - 1) * gamma(nu_val / 2)
    m_Y   <- (m_num / m_den) * (gamma - 1/gamma)

    # Second raw moment of Y (E[Y^2])
    A2 <- (nu_val / (nu_val - 2)) * (gamma^2 + 1/gamma^2 - 1)

    # Third raw moment of Y (E[Y^3])
    A3_num <- nu_val^(3/2) * gamma((nu_val - 3) / 2)
    A3_den <- sqrt(pi) * gamma(nu_val / 2)
    A3     <- (A3_num / A3_den) * ((gamma^4 - gamma^(-4)) / (gamma + 1/gamma))

    # Convert raw moments to central moments: X = Y - m_Y
    EX2 <- A2 - m_Y^2
    EX3 <- A3 - 3 * m_Y * A2 + 2 * m_Y^3

    return(EX3 / (EX2^(1.5)))
  }
}
