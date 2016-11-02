# Functions for doing connectivity estimation with multiple sites using multinomial distributions

# Function for dirichlet distribution. Adapted from the ddirichlet function in gtools
ddirichlet <- function (x, alpha, log=FALSE) 
{
  dirichlet1 <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- (alpha - 1) * log(x)
    s <- ifelse(alpha == 1 & x == 0, -Inf, s)
    return(sum(s) - logD)
  }
  if (!is.matrix(x)) 
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    else x <- t(x)
  if (!is.matrix(alpha)) 
    alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                    byrow = TRUE)
  if (any(dim(x) != dim(alpha))) 
    stop("Mismatch between dimensions of 'x' and 'alpha'.")
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
                                                         ])
  pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- -Inf
  pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- -Inf
  
  if (!log)
    pd=exp(pd)
  
  pd
}

#' Calculates unnormalized probability density for relative connectivity values 
#' from multiple distinct sites
#' 
#' This functions calculates the unnormalized probability density function for 
#' the relative (to all settlers at the destination site) connectivity value for
#' larval transport between multiple source sites to a destination site. An 
#' arbitrary number of source sites can be evaluated.
#' 
#' As this function returns the unnormalized probability density, it must be
#' normalized somehow to be produce a true probability density.  This can be
#' acheived using a variety of approaches, including brute force integration of
#' the unnormalized probability density and MCMC algorithms.
#' 
#' @param phis Vector of fractions of individuals (i.e., eggs) from the source 
#'   populations settling at the destination population
#' @param ps Vector of fractions of individuals (i.e., eggs) marked in each of 
#'   the source populations
#' @param ks Vector of numbers of marked settlers from each source population 
#'   found in the sample
#' @param n.sample Vector of total numbers of settlers collected
#' @param log Boolean indicating whether or not to return the log probability 
#'   density.  Defaults to \code{FALSE}.
#' @param dirichlet.prior.alphas Parameter value for a Dirichlet prior 
#'   distribution for the \code{phis}. Can be a single value for a Dirichlet 
#'   prior with uniform parameters, or a vector of length = 
#'   \code{length(phis)+1}. Defaults to \code{1/(length(phis)+1)}, the value for
#'   the "reference distance" non-informative prior of Berger et al. 2015.
#'   
#' @return The unnormalized probability density value.  If \code{log=TRUE}, then
#'   the logarithm of the probability density value will be returned.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @references Berger JO, Bernardo JM, Sun D (2015) Overall Objective Priors. 
#'   Bayesian Analysis 10:189-221. doi:10.1214/14-BA915
#'   
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.connectivity_estimation.multinomial.R
#' @export
d.rel.conn.multinomial.unnorm <- function(phis,ps,ks,n.sample,log=FALSE,
                                          dirichlet.prior.alphas=1/(length(phis)+1)) {
  if (sum(ks)>n.sample)
    stop("sum(ks) must be less than n")
  
  if (any(ps>1 | ps<0))
    stop("ps must be between 0 and 1")
  
  if (!all(length(phis) == c(length(ps),length(ks))))
    stop("phis, ps and ks must be vectors of same length")
  
  alphas = c(ks,n.sample-sum(ks))+1
  x = c(phis*ps,1-sum(phis*ps))
  phis = c(phis,1-sum(phis))
  
  if (length(dirichlet.prior.alphas)==1)
    dirichlet.prior.alphas=rep(dirichlet.prior.alphas,length(alphas))
  
  d = ddirichlet(x,alphas,log=T)+ddirichlet(phis,dirichlet.prior.alphas,log=T)
  if (!log)
    d = exp(d)

  return(d)
}
