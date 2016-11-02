#' Tools for Working with Connectivity Data
#' 
#' Collects several different methods for analyzing and
#' working with connectivity data in R.  Though primarily oriented
#' towards marine larval dispersal, many of the methods are general
#' and useful for terrestrial systems as well.
#'
#' \tabular{ll}{
#' Package: \tab ConnMatTools\cr
#' Type: \tab Package\cr
#' Version: \tab 0.3.3\cr
#' Date: \tab 2016-11-01\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab no\cr
#' }
#'
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @author Marco Andrello \email{marco.andrello@@gmail.com}
#' 
#' @name ConnMatTools
#' @docType package
#' @title Tools for working with connectivity matrices
#' @keywords package
#' @examples
#' \dontrun{optimalSplitConnMat(CM)}
#' @seealso See \code{\link{optimalSplitConnMat}}, \code{\link{d.rel.conn.beta.prior}}
#' 
#' @references Jacobi, M. N., and Jonsson, P. R. 2011. Optimal networks of 
#'   nature reserves can be found through eigenvalue perturbation theory of the 
#'   connectivity matrix. Ecological Applications, 21: 1861-1870.
#' @references Jacobi, M. N., Andre, C., Doos, K., and Jonsson,
#' P. R. 2012. Identification of subpopulations from connectivity
#' matrices. Ecography, 35: 1004-1016.
#' @references Gruss, A., Kaplan, D. M., and Lett, C. 2012. Estimating local 
#'   settler-recruit relationship parameters for complex spatially explicit 
#'   models. Fisheries Research, 127-128: 34-39.
#' @references Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal 
#'   per recruit: An efficient method for assessing sustainability in marine 
#'   reserve networks. Ecological Applications, 16: 2248-2263.
#' @references Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
#'   (in press) Uncertainty in empirical estimates of marine larval connectivity. 
#'   ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.
#' @references White, J. W. 2010. Adapting the steepness parameter from 
#'   stock-recruit curves for use in spatially explicit models. Fisheries 
#'   Research, 102: 330-334.
#' @references Gruss A, Kaplan DM, Hart DR (2011) Relative Impacts of Adult
#'   Movement, Larval Dispersal and Harvester Movement on the Effectiveness of
#'   Reserve Networks. PLoS ONE 6:e19960
#' @references Beverton RJH, Holt SJ (1957) On the dynamics of exploited fish
#'   populations. H.M.S.O., London. 533 pp.
#'
NULL
