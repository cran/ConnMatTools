\name{subpopsVectorToList}
\alias{subpopsVectorToList}
\title{Convert subpopulation vector to a list of indices}
\usage{
subpopsVectorToList(x)
}
\arguments{
  \item{x}{vector of subpopulation identifications}
}
\value{
A list where each element is a vector of indices for a
given subpopulation.
}
\description{
A helper function to convert a vector of subpopulation
identifications into a list appropriate for
\code{\link{recSplitConnMat}},
\code{\link{qualitySubpops}}, etc.
}
\details{
Note that subpopulations list will be ordered according to
the numerical order of the subpopulation indices in the
original matrix, which will not necessarily be that of the
spatial order of sites in the original connectivity matrix.
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\seealso{
See also \code{\link{recSplitConnMat}},
\code{\link{qualitySubpops}}
}

