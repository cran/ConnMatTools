\encoding{UTF-8}
\name{protectedAreaSelection}
\alias{protectedAreaSelection}
\title{Function to select optimal network of protected areas based on connectivity}
\usage{
protectedAreaSelection(conn.mat, nev = dim(conn.mat)[1], delta = 0.1,
  theta = 0.05, M = 20, epsilon.lambda = 1e-04, epsilon.uv = 0.05,
  only.list = T, ...)
}
\arguments{
  \item{conn.mat}{a square connectivity matrix.}

  \item{nev}{number of eigenvalues and associated
  eigenvectors to be calculated.}

  \item{delta}{the effect of protecting site i (e.g.
  increase in survival or fecundity in protected areas
  relative to unprotected areas). Now a single value, in
  future it will be possible to specify site-specific
  values. The perturbation theory used in the construction
  of the algorithm assumes delta to be small (e.g.
  delta=0.1). However, higher values give also good
  results.}

  \item{theta}{the threshold of donor times recipient value
  that a site must have to be selected.}

  \item{M}{the maximal number of sites selected from each
  subpopulation even if there are more sites above the
  threshold theta}

  \item{epsilon.lambda}{Threshold for removing complex
  eigenvalues.}

  \item{epsilon.uv}{Threshold for removing eigenvectors
  with elements of opposite signs of comparable magnitude.}

  \item{only.list}{Logical, whether the function return
  only the list of selected sites or also the predicted
  impact of each selected site on the eigenvalues}

  \item{\dots}{Additional arguments for the
  \code{\link{eigs}} function.}
}
\value{
If only.list is \code{TRUE}, just returns the list of
selected sites.  If \code{FALSE}, then result will be a
list containing selected sites and predicted impact of each
selected site on the eigenvalues.
}
\description{
This function finds the optimal network of protected areas
based on connectivity using the eigenvalue perturbation
approach described in Nilsson Jacobi & Jonsson (2011).
}
\author{
Marco Andrello \email{marco.andrello@gmail.com}
}
\references{
Jacobi, M. N., and Jonsson, P. R. 2011. Optimal networks of
nature reserves can be found through eigenvalue
perturbation theory of the connectivity matrix. Ecological
Applications, 21: 1861–1870.
}

