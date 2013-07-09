\name{localRetention}
\alias{localRetention}
\title{Local retention of a connectivity matrix}
\usage{
localRetention(conn.mat)
}
\arguments{
  \item{conn.mat}{A square connectivity matrix.}
}
\description{
Local retention is defined as the diagonal elements of the
connectivity matrix.
}
\examples{
library(ConnMatTools)
data(chile.loco)

sr <- selfRecruitment(chile.loco)
lr <- localRetention(chile.loco)
rlr <- relativeLocalRetention(chile.loco)
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}

