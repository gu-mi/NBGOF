
################################################################################
# helper functions in NBGOF, including discarded functions for future use
# NOTE: all functions should be made invisible to end-users
################################################################################

#' @title Group-wise sort of residual matrix in multivariate case
#' @description Group-wise sort of residual matrix in multivariate case.
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
sort.vec = function(x, grp.ids)  ave(x, grp.ids, FUN = sort)

printf = function(...) cat(sprintf(...));
  
#' @title Based on p-values obtained, perform chi-square-like tests for overall GOF
#' 
#' @description This function is designed to provide different methods of combining p-values
#' 
#' @param x an object of class \code{gofm}
#' @param ... for future use
#' 
#' @return P-values from different methods
#' 
#' @details When the response is a count matrix, we can use this function to test
#' the goodness-of-fit of a specified negative binomial dispersion model. 
#' 
#' @usage chisq_gof(x, ...)
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
chisq_gof = function(x, ...){
  
  ## extract quantities from "gofm" object
  sim = x$sim
  v.pvals = x$v.pvals  # p-values based on vertical distances
  p.pvals.1s = x$p.pvals.1s  # p-values based on Pearson statistics (1-sided)
  m = length(v.pvals)
  
  ## Fisher's method of combining p-values
  # vertical distance measure
  Fisher.method.X2.vert = -2*sum(log(v.pvals))
  Fisher.method.p.vert = pchisq(Fisher.method.X2.vert, 2*m, lower.tail=FALSE) 
  
  # Pearson statistics (1-sided)
  Fisher.method.X2.pear.1s = -2*sum(log(p.pvals.1s))
  Fisher.method.p.pear.1s = pchisq(Fisher.method.X2.pear.1s, 2*m, lower.tail=FALSE) 
  
  # save as a list
  results = list(fisher.vert = Fisher.method.p.vert,
                 fisher.pear.1s = Fisher.method.p.pear.1s
  )
  return(results)
  
}
