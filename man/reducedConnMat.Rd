\encoding{UTF-8}
\name{reducedConnMat}
\alias{reducedConnMat}
\title{Reduced connectivity matrix according to a set of subpopulations}
\usage{
reducedConnMat(subpops.lst, conn.mat)
}
\arguments{
  \item{subpops.lst}{A list whose elements are vectors of
  indices for each subpopulation.  If a vector of integers
  is given, then \code{\link{subpopsVectorToList}} is
  applied to convert it to a list of subpopulations.}

  \item{conn.mat}{A square connectivity matrix.}
}
\value{
A reduced connectivity matrix.  The sum of all elements of
this reduced connectivity matrix will be equal to the sum
of all elements of the original connectivity matrix.
}
\description{
Reduces a connectivity matrix based on a set of
subpopulations.  If there are N subpopulations, then the
reduced matrix will have dimensions NxN.  The reduced
matrix will be ordered according to the order of
subpopulations in \code{subpops.lst}.
}
\examples{
library(ConnMatTools)
data(chile.loco)

num <- prod(dim(chile.loco)) / sum(chile.loco)
betas <- betasVectorDefault(n=num,steps=4)
chile.loco.split <- optimalSplitConnMat(chile.loco,normalize.cols=FALSE,
                                        betas=betas)

# Extra 3rd division
print(paste("Examining split with",names(chile.loco.split$best.splits)[3],
            "subpopulations."))
pops <- subpopsVectorToList(chile.loco.split$subpops[,chile.loco.split$best.splits[[3]]$index])

reduce.loco <- reducedConnMat(pops,chile.loco)

sr <- selfRecruitment(reduce.loco)
lr <- localRetention(reduce.loco)
rlr <- relativeLocalRetention(reduce.loco)
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
See also \code{\link{qualitySubpops}}
}

