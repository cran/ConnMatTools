################################################
# Keeping this bootstrap approach around, but not exporting it
################################################

# Frequency of the number of marked settlers for finite total settlement
# estimated using a bootstrap approach
# 
# This function calculates the relative frequency with which a number of marked 
# settlers would be observed as a function of relative connectivity, sample size
# and the total number of settlers to a site. Probable frequency of observation
# is calculated by randomly drawing \code{n.bootstraps} artificial samples and
# calculating for each the number of would-be marked settlers in the sample for
# all values of \code{n.origin}, the true number of settlers from the marked
# site in the entire cohort. \code{n.origin} is related to relative 
# connectivity, \code{phi}, by \code{phi=n.origin/n.settlers}.
# 
# @param p Fraction of individuals (i.e., eggs) marked in the source population
# @param n.obs Total number of settlers collected
# @param n.settlers Total number of settlers at the destination site from which 
#   the \code{n.obs} (\code{<=n.settlers}) settlers are collected
# @param n.bootstraps Number of simulated samples to generate to calculate 
#   relative frequencies. Defaults to \code{1000}.
# @param n.origin Vector of integers of possible numbers of settlers in the 
#   cohort that originated at the site of marking. All values should be integers
#   \code{<=n.settlers}.  Defaults to \code{0:n.settlers}
#   
# @return A list with the following elements: \describe{\item{res}{matrix with 
#   \code{n.obs+1} rows and \code{length(n.origin)} columns, the element of 
#   which are counts indicating the number of times among the 
#   \code{n.bootstraps} samples that \code{k} marked settlers are found among 
#   the \code{n.obs} settlers collected from the total settlement pool of 
#   \code{n.settlers} individuals. Results for \code{k} marked settlers are in 
#   row \code{k+1} of the matrix, the first row being used for \code{k=0}.  The 
#   row names of the matrix reflect the value of \code{k}, whereas columns names
#   reflect the elements of \code{n.origin}.}\item{phi}{relative connectivity
#   values corresponding to each column of \code{res}}}
#   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#   
# @family connectivity estimation
# @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
# @encoding UTF-8
# @example tests/test.connectivity_estimation.R
# @export
# @importFrom stats runif
d.rel.conn.finite.settlement.BS <- function(p,n.obs,n.settlers,n.bootstraps=1000,
                                         n.origin=0:n.settlers) {
  res = matrix(0,nrow=n.obs+1,ncol=length(n.origin))
  rownames(res) = 0:n.obs
  colnames(res) = n.origin
  
  for (k in 1:n.bootstraps) {
    x = sort(sample(n.settlers,n.obs))
    y = as.numeric(runif(n.obs)<=p)
    i=j=0
    for (l in 1:length(n.origin)) {
      if (i < n.obs) {
        for (i in (i+1):n.obs) {
          if (x[i]>n.origin[l]) {
            i = i-1
            break            
          }
          j = j + y[i]
        }
      }
      res[j+1,l] = res[j+1,l]+1
    }
    
  }
  
  return(list(res=res,n.origin=n.origin,phi=n.origin/n.settlers))
}



################################################
# Alternative using Hypergeometric distribution to calculate probability
# Better, analytic, faster
################################################

# Helper function to make sure input arguments make sense
fs.checkparams = function(n.origin,p,k,n.obs,n.settlers,prior.n.origin,q) {
  # Check for issues with input parameters
  if (!isTRUE(all.equal(n.origin,round(n.origin))) | any(n.origin>n.settlers) | 
      any(n.origin<0))
    stop("n.origin values must be integers >=0 and <=n.settlers")
  n.origin=round(n.origin)
  
  if (!isTRUE(all.equal(k,round(k))) | any(k>n.obs) | any(k<0))
    stop("k value must be an integer >=0 and <=n.obs")
  k=round(k)
  
  if (!isTRUE(all.equal(n.obs,round(n.obs))))
    stop("n.obs value must be an integer")
  n.obs=round(n.obs)
  
  if (!isTRUE(all.equal(n.settlers,round(n.settlers))) | any(n.settlers<n.obs))
    stop("n.settlers value must be an integer >= n.obs")
  n.settlers=round(n.settlers)
  
  if (length(prior.n.origin) != 1 & length(prior.n.origin) != n.settlers+1)
    stop("prior.n.origin must be either a scalar or a vector of length n.settlers+1")
  
  if (any(p<0) | any(p>1))
    stop("Value of p must be >=0 and <=1")
  
  if (any(q<0) | any(q>1))
    stop("q values must be >=0 and <=1")
  
  return(list(n.origin=n.origin,p=p,k=k,n.obs=n.obs,n.settlers=n.settlers,
              prior.n.origin=prior.n.origin,q=q))
}

#' Estimate the probability distribution for the number of settlers originating 
#' at a site given a sample from a finite settler pool
#' 
#' These functions calculate the probability mass function 
#' (\code{d.rel.conn.finite.settlement}), the cumulative distribution function 
#' (\code{p.rel.conn.finite.settlement}) and the quantile function 
#' (\code{q.rel.conn.finite.settlement}) for the true number of settlers at a 
#' site that originated in a particular site given a known fraction of marked 
#' eggs among the eggs originating at the source site, a sample of settlers at 
#' the destination site, a known fraction of which are marked, and a finite 
#' settler pool of known size.
#' 
#' The relative connectivity between the source and destination sites is 
#' calculated as \code{n.origin/n.settlers}.
#' 
#' @param n.origin Vector of integers of possible numbers of settlers in the 
#'   cohort that originated at the site of marking. All values should be 
#'   integers \code{<=n.settlers}.
#' @param q Vector of quantiles
#' @param p Fraction of individuals (i.e., eggs) marked in the source population
#' @param k Number of marked settlers in sample
#' @param n.obs Total number of settlers collected
#' @param n.settlers Total number of settlers at the destination site from which
#'   the \code{n.obs} (\code{<=n.settlers}) settlers are collected
#' @param prior.n.origin A prior probability mass function for the number of 
#'   settlers in the cohort originating at the site of marking. Must be a scalar
#'   or a vector of length \code{n.settlers+1}. Defaults to \code{1}.
#'   
#' @return A vector of probabilities or quantiles.
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @describeIn d.rel.conn.finite.settlement Returns the probability mass 
#'   function for the numbers of settlers in the cohort that originated at the
#'   source site (i.e., site of marking).
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @example tests/test.connectivity_estimation.R
#' @export
#' @importFrom stats dbinom
#' @importFrom stats dhyper
d.rel.conn.finite.settlement <- function(n.origin,p,k,n.obs,n.settlers,
                                  prior.n.origin=1) {
  l=fs.checkparams(n.origin=n.origin,p=p,k=k,n.obs=n.obs,n.settlers=n.settlers,
                   prior.n.origin=prior.n.origin,q=0)
  for (nn in names(l))
    assign(nn,l[[nn]])
  
  # Precalculate hypergeometric stuff for speed
  dh = dhyper(k,0:n.settlers,n.settlers-(0:n.settlers),n.obs) 
  # n.origin=0 value corresponds to first entry
  
  f = function(N) {
    r = sapply(N,function(NN) sum(dbinom(k:NN,NN,p)*dh[(k:NN)+1]) )
    return(r)
  }
  
  ff = f(0:n.settlers) * prior.n.origin
  
  ff[n.origin+1] / sum(ff)
}

#' @describeIn d.rel.conn.finite.settlement Returns the cumulative distribution 
#'   function for the numbers of settlers in the cohort that originated at the
#'   source site (i.e., site of marking).
#' @export
p.rel.conn.finite.settlement <- function(n.origin,p,k,n.obs,n.settlers,
                                         prior.n.origin=1) {
  l=fs.checkparams(n.origin=n.origin,p=p,k=k,n.obs=n.obs,n.settlers=n.settlers,
                   prior.n.origin=prior.n.origin,q=0)
  for (nn in names(l))
    assign(nn,l[[nn]])

  ff = d.rel.conn.finite.settlement(0:max(n.origin),p,k,n.obs,n.settlers,prior.n.origin)
  fff = cumsum(ff)

  return(fff[n.origin+1])
}

#' @describeIn d.rel.conn.finite.settlement Returns quantiles of the cumulative
#'   distribution function for the numbers of settlers in the cohort that
#'   originated at the source site (i.e., site of marking).
#' @export
q.rel.conn.finite.settlement <- function(q,p,k,n.obs,n.settlers,
                                         prior.n.origin=1) {
  l=fs.checkparams(n.origin=0,p=p,k=k,n.obs=n.obs,n.settlers=n.settlers,
                   prior.n.origin=prior.n.origin,q=q)
  for (nn in names(l))
    assign(nn,l[[nn]])
  
  ff = p.rel.conn.finite.settlement(0:n.settlers,p,k,n.obs,n.settlers,prior.n.origin)
  
  # Special code to agree with behavior of qhyper and qbinom
  # Works, but there has to be a better way
  
  # Remove any zero probabilities from PMS
  I = which(ff>0)
  ff2=ff[I]
  
  # Now test for all values
  # Can do this without further testing because we have made sure q>=0 and q<=1
  I2 = sapply(q,function(qq) min(which(qq<=ff2)) )
  
  I1 = I[I2]-1
  return(I1)  
}
