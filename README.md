# ConnMatTools - Tools for working with connectivity matrices

## Description

Collects several different methods for analyzing and
working with connectivity matrices in R.  Though primarily oriented
towards marine larval dispersal, many of the methods are general and
useful for terrestrial systems as well.

## Maintainer

[David M. Kaplan](mailto:dmkaplan2000@gmail.com)

## Installation of latest version via github

The latest development version of the package can be installed directly from
[github](https://github.com/dmkaplan2000/ConnMatTools) using the
[devtools](https://cran.r-project.org/package=devtools)
R package. First that package needs to be installed:

     install.packages("devtools")

Once this is done, one can install the development version with:

     require("devtools")
     install_github("ConnMatTools","dmkaplan2000","master")

Replace "master" with another branch, tag or commit to obtain a
different version of the package.

## References

Methods from the following publications are implemented in the package:

Jacobi, M. N., and Jonsson, P. R. 2011. Optimal networks of 
  nature reserves can be found through eigenvalue perturbation theory of the 
  connectivity matrix. Ecological Applications, 21: 1861–1870.
  
Jacobi, M. N., André, C., Döös, K., and Jonsson,
P. R. 2012. Identification of subpopulations from connectivity
matrices. Ecography, 35: 1004-1016.

Grüss, A., Kaplan, D. M., and Lett, C. 2012. Estimating local 
  settler-recruit relationship parameters for complex spatially explicit 
  models. Fisheries Research, 127-128: 34-39.

Kaplan, D. M., Botsford, L. W., and Jorgensen, S. 2006. Dispersal 
  per recruit: An efficient method for assessing sustainability in marine 
  reserve networks. Ecological Applications, 16: 2248-2263.

Kaplan DM, Cuif M, Fauvelot C, Vigliola L, Nguyen-Huu T, Tiavouane J and Lett C 
  (in press) Uncertainty in empirical estimates of marine larval connectivity. 
  ICES Journal of Marine Science. doi:10.1093/icesjms/fsw182.

White, J. W. 2010. Adapting the steepness parameter from 
  stock-recruit curves for use in spatially explicit models. Fisheries 
  Research, 102: 330-334.

Grüss A, Kaplan DM, Hart DR (2011) Relative Impacts of Adult
  Movement, Larval Dispersal and Harvester Movement on the Effectiveness of
  Reserve Networks. PLoS ONE 6:e19960

Beverton RJH, Holt SJ (1957) On the dynamics of exploited fish
  populations. H.M.S.O., London. 533 pp.
