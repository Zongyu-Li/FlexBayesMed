% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dctpt.R
\name{dctpt}
\alias{dctpt}
\title{Density of the Centered Two-Piece Student's \eqn{t}-Distribution (CTPT)}
\usage{
dctpt(x, gamma, nu, log = FALSE)
}
\arguments{
\item{x}{Numeric vector. Quantiles at which to compute the density.}

\item{gamma}{Positive numeric scalar. Skewness parameter.}

\item{nu}{Numeric scalar \eqn{\nu > 1} (degrees of freedom) or "normal".}

\item{log}{Logical. Return log-density? Default: \code{FALSE}.}
}
\value{
Numeric vector of densities.
}
\description{
Density of the Centered Two-Piece Student's \eqn{t}-Distribution (CTPT)
}
\details{
The CTPT distribution is defined by shifting the original two-piece Student's \eqn{t}-distribution
(with skewness parameter \eqn{\gamma} and degrees of freedom \eqn{\nu}) by its mean, ensuring the
resulting distribution has zero mean. The density of the CTPT distribution is given by:
\deqn{%
  \begin{aligned}
      p(\varepsilon|\gamma,\nu) = \frac{2}{\gamma + \frac{1}{\gamma}} \Bigg[& f_\nu\left(\frac{\varepsilon + m(\gamma,\nu)}{\gamma}\right) I_{[-m(\gamma,\nu),+\infty)}(\varepsilon) \;+ \\
      &f_\nu\left(\gamma\big(\varepsilon + m(\gamma,\nu)\big)\right) I_{(-\infty,-m(\gamma,\nu))}(\varepsilon) \Bigg], \quad \gamma \in \mathbb{R}_+, \; \nu > 1,
  \end{aligned}
}
where:
\itemize{
\item \eqn{f_\nu} is the Student's \eqn{t}-density with \eqn{\nu} degrees of freedom.
\item \eqn{m(\gamma,\nu)} is the mean of the unadjusted two-piece distribution, computed by \code{\link{mean_two_piece_student}}.
\item \eqn{I_A(\varepsilon)} denotes the indicator function for set \eqn{A}.
}
}
\examples{
x <- seq(-5, 5, length.out = 100)
density <- dctpt(x, gamma = 2, nu = 5)
plot(x, density, type = "l")
}
