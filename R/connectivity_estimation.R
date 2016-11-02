################################################################
# This code estimates probability distribution for a 
# single type of marked individuals among a set of collected individuals (larvae) 
# from a source population with known fraction of marked individuals (eggs).
# 
# In the code below, I first estimate probability distributions assuming a
# uniform prior for connectivitt.  
#
# Then I do the more complicated case of a general prior.
################################################################

################################################
# Functions for calculating PDF for uniform prior
################################################

#' Estimate the probability distribution of relative connectivity values 
#' assuming a uniform prior distribution
#' 
#' These functions calculate the probability density function 
#' (\code{d.rel.conn.unif.prior}), the probability distribution function (aka 
#' the cumulative distribution function; \code{p.rel.conn.unif.prior}) and the 
#' quantile function (\code{q.rel.conn.unif.prior}) for the relative (to all 
#' settlers at the destination site) connectivity value for larval transport 
#' between a source and destination site given a known fraction of marked 
#' individuals (i.e., eggs) in the source population.  A uniform prior is used
#' for the relative connectivity value.
#' 
#' Estimations of the probability distribution are derived from the Beta 
#' distribution (see \code{\link{dbeta}}) and should be exact to great 
#' precision.
#' 
#' @param phi Vector of fractions of individuals (i.e., eggs) from the source 
#'   population settling at the destination population
#' @param q Vector of quantiles
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param k Number of marked settlers found in sample
#' @param n Total number of settlers collected
#' @param log If \code{TRUE}, returns natural logarithm of probabilities, except
#'   for \code{\link{q.rel.conn.unif.prior}}, which expects log of probabilities
#'   as inputs
#' @param \dots Extra arguments to Beta distribution functions.  See 
#'   \code{\link{dbeta}} for details.  For expert use only.
#'   
#' @return Vector of probabilities or quantiles.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @describeIn d.rel.conn.unif.prior Returns the probability density for 
#'   relative connectivity between a pair of sites
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @example tests/test.connectivity_estimation.R
#' @export
#' @importFrom stats dbeta
#' @importFrom stats pbeta
d.rel.conn.unif.prior <- function(phi,p,k,n,log=FALSE,...) {
  if (log) {
    log(p) + dbeta(phi*p,k+1,n-k+1,log=log,...) - pbeta(p,k+1,n-k+1,log.p=log,...)
  } else {
    p*dbeta(phi*p,k+1,n-k+1,log=log,...)/pbeta(p,k+1,n-k+1,log.p=log,...) 
  } 
}

#' @describeIn d.rel.conn.unif.prior Returns the cumulative probability
#'   distribution for relative connectivity between a paire of sites
#' @export
#' @importFrom stats dbeta
#' @importFrom stats pbeta
p.rel.conn.unif.prior <- function(phi,p,k,n,log=FALSE,...) {
  if (log) {
    pbeta(phi*p,k+1,n-k+1,log.p=log,...) - pbeta(p,k+1,n-k+1,log.p=log,...)
  } else {
    pbeta(phi*p,k+1,n-k+1,log.p=log,...) / pbeta(p,k+1,n-k+1,log.p=log,...)
  } 
}

#' @describeIn d.rel.conn.unif.prior Estimates quantiles for the probability
#'   distribution function for relative connectivity between a pair of sites
#' @export
#' @importFrom stats dbeta
#' @importFrom stats pbeta
#' @importFrom stats qbeta
q.rel.conn.unif.prior <-function(q,p,k,n,log=FALSE,...) {
  if (log) {
    qbeta(q + pbeta(p,k+1,n-k+1,log.p=log,...),k+1,n-k+1,log.p=log,...)/p    
  } else {
    qbeta(q*pbeta(p,k+1,n-k+1,log.p=log,...),k+1,n-k+1,log.p=log,...)/p
  }
}


################################################
# Functions for calculating PDF for Beta prior
################################################

#' Estimate the probability distribution of relative connectivity values 
#' assuming a Beta-distributed prior
#' 
#' These functions calculate the probability density function 
#' (\code{d.rel.conn.beta.prior}), the probability distribution function (aka 
#' the cumulative distribution function; \code{p.rel.conn.beta.prior}) and the 
#' quantile function (\code{q.rel.conn.beta.prior}) for the relative (to all 
#' settlers at the destination site) connectivity value for larval transport 
#' between a source and destination site given a known fraction of marked 
#' individuals (i.e., eggs) in the source population.  A non-uniform prior is 
#' used for the relative connectivity value.
#' 
#' The prior distribution for relative connectivity \code{phi} defaults to a 
#' Beta distribution with both shape parameters equal to 0.5.  This is the 
#' Reference or Jeffreys prior for a binomial distribution parameter.  Both 
#' shape parameters equal to 1 corresponds to a uniform prior.
#' 
#' Estimations of the probability distribution are based on numerical 
#' integration using the \code{\link{integrate}} function, and therefore are 
#' accurate to the level of that function.  Some modification of the default 
#' arguments to that function may be necessary to acheive good results for 
#' certain parameter values.
#' 
#' @param phi Vector of fractions of individuals (i.e., eggs) from the source 
#'   population settling at the destination population
#' @param q Vector of quantiles
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param k Number of marked settlers found in sample
#' @param n Total number of settlers collected
#' @param prior.shape1 First shape parameter for Beta distributed prior. 
#'   Defaults to 0.5.
#' @param prior.shape2 Second shape parameter for Beta distributed prior. 
#'   Defaults to being the same as \code{prior.shape1}.
#' @param prior.func Function for prior distribution.  Should take one 
#'   parameter, \code{phi}, and return a probability.  Defaults to 
#'   \code{function(phi) dbeta(phi,prior.shape1,prior.shape2)}.  If this is
#'   specified, then inputs \code{prior.shape1} and \code{prior.shape2} are
#'   ignored.
#' @param N Number of points at which to estimate cumulative probability 
#'   function for reverse approximation of quantile distribution. Defaults to 
#'   \code{1000}.
#' @param \dots Extra arguments for the \code{\link{integrate}} function used 
#'   for normalization of probability distributions.
#'   
#' @return Vector of probabilities or quantiles, or a function in the case of 
#'   \code{\link{q.rel.conn.beta.prior.func}}.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @describeIn d.rel.conn.beta.prior Returns the probability density for 
#'   relative connectivity between a pair of sites
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @example tests/test.connectivity_estimation.beta.prior.R
#' @export
#' @importFrom stats dbeta
d.rel.conn.beta.prior <- function(phi,p,k,n,
                                  prior.shape1=0.5,
                                  prior.shape2=prior.shape1,
                                  prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                  ...) {
  f = function(phi) dbeta(phi*p,k+1,n-k+1) * prior.func(phi)
  ff = integrate(f,0,1,...)$value
  
  f(phi) / ff
}

#' @describeIn d.rel.conn.beta.prior Returns the cumulative probability
#'   distribution for relative connectivity between a paire of sites
#' @export
#' @importFrom stats dbeta
#' @importFrom stats integrate
p.rel.conn.beta.prior <- function(phi,p,k,n,
                                  prior.shape1=0.5,
                                  prior.shape2=prior.shape1,
                                  prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                  ...) {
  f = function(phi) dbeta(phi*p,k+1,n-k+1) * prior.func(phi)
  ff = integrate(f,0,1,...)$value
  fff = function(phi) {
    if (phi==0) return(0)
    
    integrate(f,0,phi,...)$value
  }
  
  x = sapply(phi,fff) / ff
  
  # Cludge to fix imprecise integration leading to values greater than 1
  x[x>1 & x<(1.000001)]=1
  
  return(x)
}

#' @describeIn d.rel.conn.beta.prior Returns a function to estimate quantiles for
#'   the probability distribution function for relative connectivity between a
#'   pair of sites.
#' @include utils.R
#' @export
#' @importFrom stats dbeta
q.rel.conn.beta.prior.func <-function(p,k,n,
                                      prior.shape1=0.5,
                                      prior.shape2=prior.shape1,
                                      prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                      N=1000,
                                      ...) {
  phi = seq(0,1,length.out=N)
  q = p.rel.conn.beta.prior(phi,p,k,n,prior.func=prior.func,...)
  return(rel.conn.approxfun(q,phi))
}

#' @describeIn d.rel.conn.beta.prior Estimates quantiles for the probability
#'   distribution function for relative connectivity between a pair of sites
#' @export
#' @importFrom stats dbeta
q.rel.conn.beta.prior <-function(q,p,k,n,
                                 prior.shape1=0.5,
                                 prior.shape2=prior.shape1,
                                 prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                 N=1000,
                                 ...)
  (q.rel.conn.beta.prior.func(p,k,n,prior.func=prior.func,N=N,...))(q)
