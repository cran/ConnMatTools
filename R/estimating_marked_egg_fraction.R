#' Fraction of eggs marked for male and female mark transmission
#' 
#' Estimates the fraction of eggs produced at the source site that are the 
#' result of crossing parents, one or both of which have been genotyped. Based 
#' on the assumption that probability of breeding between pairs of individuals 
#' is completely independent of whether or not one or more of those individuals 
#' was genotyped.
#' 
#' @param p.female Fraction of all adult females genotyped in the source 
#'   population
#' @param p.male  Fraction of all adult males genotyped in the source 
#'   population. Defaults to be equal to \code{p.female}
#'   
#' @return A list with the following elements: \describe{\item{prob.matrix}{2x2 
#'   matrix with probabilities for producing offspring with male or female known
#'   or unknown parents}\item{p}{fraction of all eggs produced at source site 
#'   that will come from at least one genotyped 
#'   parent}\item{p.female.known}{Fraction of eggs with a single known female 
#'   parent among all eggs that have one or more known 
#'   parents}\item{p.male.known}{Fraction of eggs with a single known male 
#'   parent among all eggs that have one or more known 
#'   parents}\item{p.two.known.parents}{Fraction of eggs with two known parents
#'   among all eggs that have one or more known parents}}
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
dual.mark.transmission <- function(p.female,p.male=p.female) {
  if (p.female>1 || p.female<0 || p.male<0 || p.male>1 )
    stop('Probabilities must lie between 0 and 1.')
  
  pf = c(1-p.female,p.female)
  pm = c(1-p.male,p.male)
  
  pp = pf %*% t(pm)
  
  rownames(pp) = paste("female",c("unknown","known"),sep=".")
  colnames(pp) = paste("male",c("unknown","known"),sep=".")
  
  p = 1 - pp["female.unknown","male.unknown"]
  
  p.male.known = pp["female.unknown","male.known"] / p
  p.female.known = pp["female.known","male.unknown"] / p
  p.two.known.parents = pp["female.known","male.known"] / p
  
  return(list(prob.matrix=pp,p=p,p.male.known=p.male.known,
              p.female.known=p.female.known,
              p.two.known.parents=p.two.known.parents))
}

#' Estimates of fraction of eggs marked accounting for variability in 
#' reproductive output
#' 
#' This function estimates the fraction of eggs "marked" at a site (where the 
#' "mark" could be micro-chemical or genetic) taking into account uncertainty in
#' female (and potentially male in the case of dual genetic mark transmission) 
#' reproductive output. It generates a set of potential values for the fraction 
#' of eggs marked assuming that reproductive output of each marked or unmarked 
#' mature individual is given by a random variable drawn from a single 
#' probability distribution with known mean and standard deviation (or 
#' equivalently coefficient of variation) \strong{and} that the numbers 
#' of marked and unmarked individuals are large enough that the central limit 
#' theorem applies and, therefore, their collective reproductive outputs are 
#' reasonably well described by a gamma distribution whose mean and standard 
#' deviation are appropriately scaled based on the number of individual 
#' reproducers.  The function also returns the total egg production 
#' corresponding to each fraction of marked eggs, needed for estimating absolute
#' connectivity values (i.e., elements of the connectivity matrix needed for 
#' assessing population persistence).
#' 
#' @param n Number of random values to estimates
#' @param n.females Total number of mature females in the population
#' @param n.marked.females Number of marked females in population
#' @param mean.female Mean egg production of each mature female. Defaults to 1.
#' @param cv.female Coefficient of variation of reproductive output of an 
#'   individual mature female
#' @param dual Logical variable. If \code{TRUE}, then the fraction of marked 
#'   eggs is calculated assuming dual (male and female) mark transmission. 
#'   Defaults to \code{FALSE}.
#' @param male.uncert Logical variable. If \code{TRUE}, then variability in male
#'   sperm output is also taken into account when estimating the number of 
#'   marked eggs. Defaults to \code{FALSE}.
#' @param n.males Total number of mature males in the population. Only used if 
#'   \code{dual=TRUE}. Defaults to being equal to \code{n.females}.
#' @param n.marked.males Number of marked males in population. Only used if 
#'   \code{dual=TRUE}. Defaults to being equal to \code{n.marked.females}.
#' @param mean.male Mean sperm production of each mature male. Only used if 
#'   \code{dual=TRUE} and \code{male.uncert=TRUE}. Defaults to being equal to 
#'   \code{mean.female}.
#' @param cv.male Coefficient of variation of reproductive output of an 
#'   individual mature male. Only used if \code{dual=TRUE} and 
#'   \code{male.uncert=TRUE}. Defaults to being equal to \code{cv.female}.
#' @param p.marked.females Fraction of marked females in population. Can be 
#'   supplied instead of \code{n.marked.females}. Ignored if 
#'   \code{n.marked.females} is given.
#' @param p.marked.males Fraction of marked males in population. Can be supplied
#'   instead of \code{n.marked.males}. Only used if \code{dual=TRUE}. Ignored if
#'   \code{n.marked.males} is given.
#'   
#' @return A list with the following elements: \describe{\item{p}{Vector of 
#'   length \code{n} with estimates for fraction of marked 
#'   eggs}\item{eggs}{Vector of length \code{n} with estimates for total egg 
#'   production}\item{marked.eggs}{Vector of length \code{n} with estimates for 
#'   total number of marked eggs produced}\item{sperm}{Only returned if 
#'   \code{dual=TRUE}. If \code{male.uncert=FALSE}, then a scalar equal to 
#'   \code{n.males}. Otherwise, a vector of length \code{n} with estimates for 
#'   total sperm production}\item{marked.sperm}{Only returned if 
#'   \code{dual=TRUE}. If \code{male.uncert=FALSE}, then a scalar equal to 
#'   \code{n.marked.males}. Otherwise, a vector of length \code{n} with
#'   estimates for total marked sperm production}}
#'   
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#'   
#' @family connectivity estimation
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @example tests/test.marked.egg.fraction.R
#' @encoding UTF-8
#' @export
#' @include utils.R
#' @importFrom stats rgamma
r.marked.egg.fraction <- function(n,
                                  n.females,
                                  n.marked.females=round(n.females*p.marked.females),
                                  mean.female=1,cv.female,
                                  dual=FALSE,male.uncert=FALSE,
                                  n.males=n.females,
                                  n.marked.males=tryCatch(round(n.males*p.marked.males),error=function(e) n.marked.females),
                                  mean.male=mean.female,cv.male=cv.female,
                                  p.marked.females,p.marked.males=p.marked.females
                                  ) {
  if (n.marked.females>n.females)
    stop("n.females must be greater than n.marked.females")
  
  f1 = function(x,y) x / (x+y)
  
  lf = gammaParamsConvert(mean=mean.female,sd=mean.female*cv.female)

  xf = rgamma(n,n.marked.females*lf$shape,lf$scale)
  yf = rgamma(n,(n.females-n.marked.females)*lf$shape,lf$scale)
  
  if (!dual)
    return(list(eggs=xf+yf,marked.eggs=xf,p=f1(xf,yf)))
  
  if (n.marked.males>n.males)
    stop("n.males must be greater than n.marked.males")

  xm = n.marked.males
  ym = n.males - n.marked.males
  
  if (male.uncert) {
    lm = gammaParamsConvert(mean=mean.male,sd=mean.male*cv.male)
    xm = rgamma(n,n.marked.males*lm$shape,lm$scale)
    ym = rgamma(n,(n.males-n.marked.males)*lm$shape,lm$scale)    
  }
  
  # Function that gives fraction marked when there is dual mark transmission
  f2 = function(pf,pm) pm + pf - pm*pf
  
  return(list(eggs=xf+yf,marked.eggs=xf,
              sperm=xm+ym,marked.sperm=xm,
              p=f2(f1(xf,yf),f1(xm,ym))))
}
