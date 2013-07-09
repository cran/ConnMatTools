\name{hockeyStick}
\alias{hockeyStick}
\title{Hockey-stick settler-recruit relationship}
\usage{
hockeyStick(S, slope = 1/0.35, Rmax = 1)
}
\arguments{
  \item{S}{a vector of settlement values, 1 for each site.}

  \item{slope}{slope at the origin of the settler-recruit
  relationship.  Can be a vector of same length as
  \code{S}.}

  \item{Rmax}{maximum recruitment value.}
}
\value{
A vector of recruitment values.
}
\description{
Calculates recruitment based on a settler-recruit
relationship that increases linearly until it reaches a
maximum values.
}
\details{
\code{slope} and \code{Rmax} can both either be scalars or
vectors of the same length as \code{S}.
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\references{
Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006.
Dispersal per recruit: An efficient method for assessing
sustainability in marine reserve networks. Ecological
Applications, 16: 2248-2263.
}

