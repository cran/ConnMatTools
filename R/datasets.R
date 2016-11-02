#' Connectivity matrix for loco (Concholepas concholepas) from Chile
#' 
#' Sample connectivity matrix representing potential larval dispersal of loco 
#' (Concholepas concholepas) from Chile.  The matrix is for 89 sites along the 
#' coast of Chile and is derived from a theoretical larval transport model.
#' 
#' @format A square 89x89 matrix with real, positive elements.
#' @references Garavelli L, Kaplan DM, Colas F, Stotz W, Yannicelli B, Lett C
#'   (2014) Identifying appropriate spatial scales for marine conservation and
#'   management using a larval dispersal model: The case of Concholepas
#'   concholepas (loco) in Chile. Progress in Oceanography 124:42-53.
#'   \href{http://dx.doi.org/10.1016/j.pocean.2014.03.011}{doi:10.1016/j.pocean.2014.03.011}
#' @name chile.loco
#' @docType data
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @encoding UTF-8
#' @keywords data
NULL

#' Sample LOD score data for simulated and real parent-child pairs
#' 
#' This dataset contains both simulated and real 'log of the odds ratio' (LOD) 
#' scores for potential parent-child pairs of humbug damselfish (\emph{Dascyllus
#' aruanus}) from New Caledonia. Data was generated using 
#' \href{http://www.pierroton.inra.fr/genetics/labo/Software/Famoz/index.html}{FaMoz}.
#' In all cases, results are for the potential parent with the highest LOD score
#' for a given larval fish (child). Simulated data is based on artificial 
#' children generated from either real potential parent-pairs (the 'in' group) 
#' or artificial parents generated from observed allelic frequencies (the 'out' 
#' group).
#' 
#' @format A list with 3 elements: \describe{\item{in.group}{5000 maximum LOD 
#'   scores for simulated children from random crossing of real potential 
#'   parents}\item{out.group}{5000 maximum LOD scores for simulated children 
#'   from random crossing of artificial potential parents based on observed 
#'   allelic frequencies}\item{real.children}{Maximum LOD scores for 200 real
#'   juvenile fish}}
#' @references Gerber S, Chabrier P, Kremer A (2003) FAMOZ: a software for 
#'   parentage analysis using dominant, codominant and uniparentally inherited 
#'   markers. Molecular Ecology Notes 3:479-481. 
#'   \href{http://dx.doi.org/10.1046/j.1471-8286.2003.00439.x}{doi:10.1046/j.1471-8286.2003.00439.x}
#' @references Kaplan et al. (submitted) Uncertainty in marine larval 
#'   connectivity estimation
#' @seealso See also \code{\link{d.rel.conn.dists.func}}
#' @name damselfish.lods
#' @docType data
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @keywords data
#' @encoding UTF-8
NULL
