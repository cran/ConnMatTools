% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connectivity_estimation.R
\encoding{UTF-8}
\name{d.rel.conn.beta.prior}
\alias{d.rel.conn.beta.prior}
\alias{p.rel.conn.beta.prior}
\alias{q.rel.conn.beta.prior.func}
\alias{q.rel.conn.beta.prior}
\title{Estimate the probability distribution of relative connectivity values 
assuming a Beta-distributed prior}
\usage{
d.rel.conn.beta.prior(
  phi,
  p,
  k,
  n,
  prior.shape1 = 0.5,
  prior.shape2 = prior.shape1,
  prior.func = function(phi) dbeta(phi, prior.shape1, prior.shape2),
  ...
)

p.rel.conn.beta.prior(
  phi,
  p,
  k,
  n,
  prior.shape1 = 0.5,
  prior.shape2 = prior.shape1,
  prior.func = function(phi) dbeta(phi, prior.shape1, prior.shape2),
  ...
)

q.rel.conn.beta.prior.func(
  p,
  k,
  n,
  prior.shape1 = 0.5,
  prior.shape2 = prior.shape1,
  prior.func = function(phi) dbeta(phi, prior.shape1, prior.shape2),
  N = 1000,
  ...
)

q.rel.conn.beta.prior(
  q,
  p,
  k,
  n,
  prior.shape1 = 0.5,
  prior.shape2 = prior.shape1,
  prior.func = function(phi) dbeta(phi, prior.shape1, prior.shape2),
  N = 1000,
  ...
)
}
\arguments{
\item{phi}{Vector of fractions of individuals (i.e., eggs) from the source 
population settling at the destination population}

\item{p}{Fraction of individuals (i.e., eggs) marked in the source population}

\item{k}{Number of marked settlers found in sample}

\item{n}{Total number of settlers collected}

\item{prior.shape1}{First shape parameter for Beta distributed prior. 
Defaults to 0.5.}

\item{prior.shape2}{Second shape parameter for Beta distributed prior. 
Defaults to being the same as \code{prior.shape1}.}

\item{prior.func}{Function for prior distribution.  Should take one 
parameter, \code{phi}, and return a probability.  Defaults to 
\code{function(phi) dbeta(phi,prior.shape1,prior.shape2)}.  If this is
specified, then inputs \code{prior.shape1} and \code{prior.shape2} are
ignored.}

\item{\dots}{Extra arguments for the \code{\link{integrate}} function used 
for normalization of probability distributions.}

\item{N}{Number of points at which to estimate cumulative probability 
function for reverse approximation of quantile distribution. Defaults to 
\code{1000}.}

\item{q}{Vector of quantiles}
}
\value{
Vector of probabilities or quantiles, or a function in the case of 
  \code{\link{q.rel.conn.beta.prior.func}}.
}
\description{
These functions calculate the probability density function 
(\code{d.rel.conn.beta.prior}), the probability distribution function (aka 
the cumulative distribution function; \code{p.rel.conn.beta.prior}) and the 
quantile function (\code{q.rel.conn.beta.prior}) for the relative (to all 
settlers at the destination site) connectivity value for larval transport 
between a source and destination site given a known fraction of marked 
individuals (i.e., eggs) in the source population.  A non-uniform prior is 
used for the relative connectivity value.
}
\details{
The prior distribution for relative connectivity \code{phi} defaults to a 
Beta distribution with both shape parameters equal to 0.5.  This is the 
Reference or Jeffreys prior for a binomial distribution parameter.  Both 
shape parameters equal to 1 corresponds to a uniform prior.

Estimations of the probability distribution are based on numerical 
integration using the \code{\link{integrate}} function, and therefore are 
accurate to the level of that function.  Some modification of the default 
arguments to that function may be necessary to acheive good results for 
certain parameter values.
}
\section{Functions}{
\itemize{
\item \code{d.rel.conn.beta.prior}: Returns the probability density for 
relative connectivity between a pair of sites

\item \code{p.rel.conn.beta.prior}: Returns the cumulative probability
distribution for relative connectivity between a paire of sites

\item \code{q.rel.conn.beta.prior.func}: Returns a function to estimate quantiles for
the probability distribution function for relative connectivity between a
pair of sites.

\item \code{q.rel.conn.beta.prior}: Estimates quantiles for the probability
distribution function for relative connectivity between a pair of sites
}}

\examples{
library(ConnMatTools)

k <- 10 # Number of marked settlers among sample
n.obs <- 87 # Number of settlers in sample

p <- 0.4 # Fraction of eggs that was marked
phi <- seq(0.001,1-0.001,length.out=101) # Values for relative connectivity

# Probability distribution assuming infinite settler pool and uniform prior
drc <- d.rel.conn.unif.prior(phi,p,k,n.obs)
qrc <- q.rel.conn.unif.prior(c(0.025,0.975),p,k,n.obs) # 95\% confidence interval

# Probability distribution assuming infinite settler pool and using reference/Jeffreys prior
drp <- d.rel.conn.beta.prior(phi,p,k,n.obs)
prp <- p.rel.conn.beta.prior(phi,p,k,n.obs)
qrp <- q.rel.conn.beta.prior(c(0.025,0.975),p,k,n.obs) # 95\% confidence interval

# Make a plot of different distributions
# black = Jeffreys prior; red = uniform prior
# Jeffreys prior draws distribution slightly towards zero
plot(phi,drp,type="l",main="Probability of relative connectivity values",
     xlab=expression(phi),ylab="Probability density")
lines(phi,drc,col="red")
abline(v=qrp,col="black",lty="dashed")
abline(v=qrc,col="red",lty="dashed")
}
\references{
Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
  (in press) Uncertainty in empirical estimates of marine larval connectivity. 
  ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
}
\seealso{
Other connectivity estimation: 
\code{\link{d.rel.conn.dists.func}()},
\code{\link{d.rel.conn.finite.settlement}()},
\code{\link{d.rel.conn.multinomial.unnorm}()},
\code{\link{d.rel.conn.multiple}()},
\code{\link{d.rel.conn.unif.prior}()},
\code{\link{dual.mark.transmission}()},
\code{\link{optim.rel.conn.dists}()},
\code{\link{r.marked.egg.fraction}()}
}
\author{
David M. Kaplan \email{dmkaplan2000@gmail.com}
}
\concept{connectivity estimation}
