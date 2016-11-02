# Function used to fix problem with quantile estimation when you have many exact 0's and
# 1's at beginning and end, respectively or probabilty distribution
#' @importFrom stats approxfun
rel.conn.approxfun <- function(x,y,xmin=0,xmax=1,...) {
  if (any(x>xmax) || any(x<xmin))
    stop(paste("Values outside of supposed bounds: [",xmin,",",xmax,"]"))
  
  I = which(x > xmin & x < xmax)
  
  I1 = max(1,min(I)-1)
  I2 = min(length(x),max(I)+1)
  
  return(approxfun(x[I1:I2],y[I1:I2],...))
}

#' Gamma distribution shape and scale parameters from mean and standard 
#' deviation, or vice-versa
#' 
#' Calculates shape and scale parameters for a gamma distribution from the mean 
#' and standard deviation of the distribution, or vice-versa.  One supplies 
#' either \code{mean} and \code{sd} or \code{shape} and \code{scale} and the 
#' function returns a list with all four parameter values.
#' 
# @param mean Mean of the gamma distribution
# @param sd Standard deviation of the gamma distribution
# @param shape Shape parameter of the gamma distribution
# @param scale Scale parameter of the gamma distribution
#' @param \dots This function can be run either supplying \code{mean} and
#'   \code{sd}, or supplying \code{shape} and \code{scale}, but not both pairs
#'   of parameters.
#'   
#' @return A list with \code{mean}, \code{sd}, \code{shape} and \code{scale} 
#'   parameters of the corresponding gamma distribution.
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @export
#' @examples
#' library(ConnMatTools)
#' mn <- 1
#' sd <- 0.4
#' l <- gammaParamsConvert(mean=mn,sd=sd)
#' x <- seq(0,2,length.out=50)
#' 
#' # Plot gamma and normal distributions - for sd << mean, the two should be very close
#' plot(x,dgamma(x,l$shape,scale=l$scale),
#'      main="Normal versus Gamma distributions",type="l")
#' lines(x,dnorm(x,l$mean,l$sd),col="red")
gammaParamsConvert <- function(...) {
  l = list(...)
  
  if (length(l) != 2)
    stop("Number of input arguments should be two.")
  
  if (is.null(l$mean))
    l$mean = l$shape * l$scale
  
  if (is.null(l$sd))
    l$sd = sqrt(l$shape * l$scale^2)
  
  if (is.null(l$shape))
    l$shape = (l$mean / l$sd)^2
  
  if (is.null(l$scale))
    l$scale = l$sd^2 / l$mean
  
  return(l)
}
