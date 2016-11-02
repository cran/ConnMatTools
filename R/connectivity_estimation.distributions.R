# Code for detecting the fraction of recaptured individuals that are marked when we
# know the distributions of marked and unmarked individuals

#' Create a probability density step function from a histogram object
#' 
#' This function creates a step function from the bars in a \code{histogram} 
#' object. By default, the step function will be normalized so that it 
#' integrates to 1.
#' 
#' @param h an object of type \code{histogram}
#' @param \dots Additional arguments for the default \code{\link{stepfun}} 
#'   function.
#' @param normalize Boolean indicating whether or not to normalize the output 
#'   stepfun so that it integrates to 1. Defaults to \code{TRUE}.  If 
#'   \code{FALSE}, then the function will integrate to \code{sum(h$counts)}
#'   
#' @return A function of class \code{stepfun}.  The height of the steps will be
#'   divided by the distance between breaks and possibly the total count.
#'   
#' @seealso See also \code{\link{d.rel.conn.dists.func}},
#'   \code{\link{optim.rel.conn.dists}}.
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.connectivity_estimation.distributions.R
#' @encoding UTF-8
#' @export
#' @importFrom stats stepfun
stepfun.hist <- function(h,...,normalize=TRUE) {
  n = 1
  if (normalize)
    n = sum(h$counts)
  
  return(stepfun(x=h$breaks,y=c(0,h$counts / diff(h$breaks) / n,0),...))
}

#' Returns probability density function (PDF) for a mix of marked and unmarked 
#' individuals
#' 
#' This function returns a probability density function (PDF) for scores for a 
#' mix of marked and unmarked individuals with known fraction of marked 
#' individuals. The distributions for marked individuals and for unmarked 
#' individuals must be known.
#' 
#' @param d.unmarked A function representing the PDF of unmarked individuals. 
#'   Must be normalized so that it integrates to 1 for the function to work 
#'   properly.
#' @param d.marked A function representing the PDF of marked individuals.  Must 
#'   be normalized so that it integrates to 1 for the function to work properly.
#'   
#' @return A function representing the PDF of observations drawn from the mixed 
#'   distribution of marked and unmarked individuals.  The function takes two 
#'   arguments: \code{p.marked}, the fraction of marked individuals in the
#'   distribution; and \code{obs}, a vector of observed score values.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @seealso See also \code{\link{d.rel.conn.dists.func}},
#'   \code{\link{optim.rel.conn.dists}}.
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
d.mix.dists.func <- function(d.unmarked,d.marked) {
  return(function(p.marked,obs) (1-p.marked)*d.unmarked(obs)+
           p.marked*d.marked(obs))
}

#' Returns probability a set of observations correspond to marked individuals
#' 
#' This function returns the probability each of a set of observations 
#' corresponds to a marked individual given the distribution of scores for 
#' unmarked and marked individuals and the fraction of individuals that are 
#' marked.
#' 
#' @param obs A vector of score values for a random sample of (marked and 
#'   unmarked) individuals from the population
#' @inheritParams d.mix.dists.func
#' @param phi The fraction of settlers at the destination population that 
#'   originated at the source population. Defaults to 0.5, which would
#'   correspond to an even sample of marked and unmarked individuals.
#' @param p Fraction of individuals (i.e., eggs) marked in the source 
#'   population. Defaults to 1.
#'   
#' @return A vector of the same size as \code{obs} containing the probability 
#'   that each individual is marked
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @seealso See also \code{\link{d.rel.conn.dists.func}},
#'   \code{\link{optim.rel.conn.dists}}.
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
prob.marked <- function(obs,d.unmarked,d.marked,phi=0.5,p=1) {
  p.marked = phi * p
  p.marked*d.marked(obs) / ((1-p.marked)*d.unmarked(obs)+p.marked*d.marked(obs))
}

# Logarithm of probability for a set of observed score values given a PDF
# 
# This function returns log-probability for a set of observed score values 
# given a PDF for the distribution of scores.
# 
# @param p Parameter in the PDF distribution.  This should be a single scalar
#   value, presumably the fraction of marked individuals in the population, but
#   the function is generic.
# @param obs A vector of score values
# @param dfunc A function representing the PDF of scores.  Should be normalized
#   so that it integrates to 1 for the function to return a true 
#   log-probability.
#   
# @return The log-probability.
#   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#   
# @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
# @encoding UTF-8
log.prob <- function(p,obs,dfunc) {
  v=dfunc(p,obs)
  
  if (any(v==0))
    return(-Inf)
  
  sum(log(v))  
}

#' Maximum-likelihood estimate for relative connectivity given two distributions
#' for scores for marked and unmarked individuals
#' 
#' This function calculates the value for relative connectivity that best fits a
#' set of observed score values, a pair of distributions for marked and unmarked
#' individuals and an estimate of the fraction of eggs marked in the source 
#' population, \code{p}.
#' 
#' @param obs Vector of observed score values for potentially marked individuals
#' @inheritParams d.mix.dists.func
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param phi0 Initial value for \eqn{\phi}, the fraction of settlers at the 
#'   destination population that originated at the source population, for 
#'   \code{\link{optim}} function. Defaults to 0.5.
#' @param method Method variable for \code{\link{optim}} function. Defaults to 
#'   \code{"Brent"}.
#' @param lower Lower limit for search for fraction of marked individuals. 
#'   Defaults to 0.
#' @param upper Upper limit for search for fraction of marked individuals. 
#'   Defaults to 1.
#' @param \dots Additional arguments for the \code{\link{optim}} function.
#'   
#' @return A list with results of optimization. Optimal fraction of marked 
#'   individuals is in \code{phi} field. Negative log-likelihood is in the 
#'   \code{neg.log.prob} field. See \code{\link{optim}} for other elements of
#'   list.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.connectivity_estimation.distributions.R
#' @encoding UTF-8
#' @export
#' @importFrom stats optim
optim.rel.conn.dists <- function(obs,d.unmarked,d.marked,p=1,
                                 phi0=0.5,method="Brent",lower=0,upper=1,
                                 ...) {
  f = function(phi)
    -log.prob(p*phi,obs,d.mix.dists.func(d.unmarked,d.marked))
  
  o=optim(par=phi0,fn=f,method=method,lower=lower,upper=upper,...)
  o$phi = o$par
  o$neg.log.prob = o$value
  
  return(o)
}

# Function that returns a function that calculates likelihood after removing the
# part that corresponds to optimal value.  This helps with numerics of 
# normalizing distribution
.d.rel.conn.dists.internal <- function(obs,d.unmarked,d.marked,p=1) {
  f = d.mix.dists.func(d.unmarked,d.marked)
  ff = function(obs,phi,p)
    log.prob(phi*p,obs,f)
  
  f0 = optim.rel.conn.dists(obs,d.unmarked,d.marked,p=p)$neg.log.prob
  
  return(function(phi) exp(ff(obs,phi,p)+f0))
}

# # Function to define number of steps to use for normalization of distributions
# .d.rel.conn.dists.N <- function(obs)
#   max(100,min(5000,2*length(obs)))

#' Functions for estimating the probability distribution for relative 
#' connectivity given a pair of distributions for scores for marked and unmarked
#' individuals
#' 
#' These functions return functions that calculate the probability density 
#' function (\code{d.rel.conn.dists.func}), the probability distribution 
#' function (aka the cumulative distribution function; 
#' \code{p.rel.conn.dists.func}) and the quantile function 
#' (\code{q.rel.conn.dists.func}) for relative connectivity given a set of 
#' observed score values, distributions for unmarked and marked individuals, and
#' an estimate of the fraction of all eggs marked at the source site, \code{p}.
#' 
#' The normalization of the probability distribution is carried out using a 
#' simple, fixed-step trapezoidal integration scheme.  By default, the number of
#' steps between relative connectivity values of 0 and 1 defaults to 
#' \code{2*length(obs)} so long as that number is comprised between \code{100} 
#' and \code{5000}.
#' 
#' @param obs Vector of observed score values for potentially marked individuals
#' @inheritParams d.mix.dists.func
#' @param p Fraction of individuals (i.e., eggs) marked in the source 
#'   population. Defaults to 1.
#' @param N number of steps between 0 and 1 at which to approximate likelihood 
#'   function as input to \code{\link{approxfun}}. Defaults to
#'   \code{2*length(obs)} so long as that number is comprised between \code{100}
#'   and \code{5000}.
#' @param prior.shape1 First shape parameter for Beta distributed prior. 
#'   Defaults to 0.5.
#' @param prior.shape2 Second shape parameter for Beta distributed prior. 
#'   Defaults to being the same as \code{prior.shape1}.
#' @param prior.func Function for prior distribution.  Should take one 
#'   parameter, \code{phi}, and return a probability.  Defaults to 
#'   \code{function(phi) dbeta(phi,prior.shape1,prior.shape2)}.  If this is 
#'   specified, then inputs \code{prior.shape1} and \code{prior.shape2} are 
#'   ignored.
#' @param \dots Additional arguments for the \code{\link{integrate}} function.
#'   
#' @return A function that takes one argument (the relative connectivity for 
#'   \code{d.rel.conn.dists.func} and \code{p.rel.conn.dists.func}; the quantile
#'   for \code{q.rel.conn.dists.func}) and returns the probability density, 
#'   cumulative probability or score value, respectively. The returned function 
#'   accepts both vector and scalar input values.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @describeIn d.rel.conn.dists.func Returns a function that is PDF for relative
#'   connectivity
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.connectivity_estimation.distributions.R
#' @encoding UTF-8
#' @export
#' @importFrom stats integrate
#' @importFrom stats approxfun
#' @importFrom stats dbeta
d.rel.conn.dists.func <- function(obs,d.unmarked,d.marked,p=1,
                                  N=max(100,min(5000,2*length(obs))),
                                  prior.shape1=0.5,
                                  prior.shape2=prior.shape1,
                                  prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                  ...) {
  phi=seq(0,1,length.out=N+1)
  
  f = sapply(phi,.d.rel.conn.dists.internal(obs,d.unmarked,d.marked,p=p))
  ff = approxfun(phi,f)
  fff = function(phi) ff(phi) * prior.func(phi)
  fs = integrate(fff,0,1,...)$value
  
  return(function(phi) fff(phi) / fs)
}

#' @describeIn d.rel.conn.dists.func Returns a function that is cumulative
#'   probability distribution for relative connectivity
#' @export
#' @importFrom stats integrate
#' @importFrom stats dbeta
p.rel.conn.dists.func <- function(obs,d.unmarked,d.marked,p=1,
                                  N=max(100,min(5000,2*length(obs))),
                                  prior.shape1=0.5,
                                  prior.shape2=prior.shape1,
                                  prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                  ...) {
  f = d.rel.conn.dists.func(obs=obs,d.unmarked=d.unmarked,d.marked=d.marked,p=p,
                            N=N,prior.func=prior.func,...)
  
  ff =  function(phi) {
    if (phi==0) return(0)
    
    integrate(f,0,phi,...)$value
  }
  
  # Not efficient to use sapply and reintegrate multiple times the same thing,
  # but don't see better option
  return(function(phi) sapply(phi,ff))
}

#' @describeIn d.rel.conn.dists.func Returns a function that is quantile
#'   function for relative connectivity
#' @include utils.R
#' @export
#' @importFrom stats integrate
#' @importFrom stats dbeta
q.rel.conn.dists.func <- function(obs,d.unmarked,d.marked,p=1,
                                  N=max(100,min(5000,2*length(obs))),
                                  prior.shape1=0.5,
                                  prior.shape2=prior.shape1,
                                  prior.func=function(phi) dbeta(phi,prior.shape1,prior.shape2),
                                  ...) {
  f = p.rel.conn.dists.func(obs=obs,d.unmarked=d.unmarked,d.marked=d.marked,p=p,
                            N=N,prior.func=prior.func,...)

  phi=seq(0,1,length.out=N+1)
  
  ff = f(phi)
  ff = ff / max(ff) # Need to renormalize due to slight integration errors
  
  return(rel.conn.approxfun(ff,phi,...))
}

