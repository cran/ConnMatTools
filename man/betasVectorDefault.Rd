\encoding{UTF-8}
\name{betasVectorDefault}
\alias{betasVectorDefault}
\title{Compute vector of beta values}
\usage{
betasVectorDefault(n, steps = 10, cycles = 3/4, coeff = 0.8, pwr = 3)
}
\arguments{
  \item{n}{numerator of formula from Jacobi et al. (2012).
  Normally will be the number of columns in the
  connectivity matrix if one normalizes the columns
  (otherwise, it would typically be \code{N^2 /
  sum(conn.mat)}, where \code{N} is the number of columns
  of \code{conn.mat}.}

  \item{steps}{number of beta values to return.  Defaults
  to 10.}

  \item{cycles}{how many cycles of \code{2*pi} to do.}

  \item{coeff}{coefficient in front of sine function}

  \item{pwr}{exponent in denominator}
}
\value{
vector of beta values
}
\description{
Helper function to compute a set of beta values using
formula used in Jacobi et al. (2012).
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\references{
Jacobi, M. N., André, C., Döös, K., and Jonsson, P. R.
2012.  Identification of subpopulations from connectivity
matrices. Ecography, 35: 1004-1016.
}
\seealso{
See also \code{\link{optimalSplitConnMat}}
}

