
################################################################################
# helper functions in NBGOF, including discarded functions for future use
# NOTE: all functions should be made invisible to end-users!
################################################################################

#' @title Group-wise sort of residual matrix in multivariate case
#' @description Group-wise sort of residual matrix in multivariate case. THIS MAY BE USEFUL
#' WHEN WE CONSIDER MULTIPLE GROUPS AND COMPLEX DESIGNS: HOW WE SORT THE RESIDUALS? CURRENTLY
#' NOT USED.
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
  p.pvals.2s = x$p.pvals.2s  # p-values based on Pearson statistics (2-sided)
  m = length(v.pvals)
  
  ## Fisher's method of combining p-values
  # vertical distance measure
  Fisher.method.X2.vert = -2*sum(log(v.pvals))
  Fisher.method.p.vert = pchisq(Fisher.method.X2.vert, 2*m, lower.tail=FALSE) 
  # Pearson statistics (1-sided)
  Fisher.method.X2.pear.1s = -2*sum(log(p.pvals.1s))
  Fisher.method.p.pear.1s = pchisq(Fisher.method.X2.pear.1s, 2*m, lower.tail=FALSE) 
  # Pearson statistics (2-sided)
  Fisher.method.X2.pear.2s = -2*sum(log(p.pvals.2s))
  Fisher.method.p.pear.2s = pchisq(Fisher.method.X2.pear.2s, 2*m, lower.tail=FALSE) 
  
  # save as a list
  results = list(fisher.vert = Fisher.method.p.vert,
                 fisher.pear.1s = Fisher.method.p.pear.1s,
                 fisher.pear.2s = Fisher.method.p.pear.2s
                 )
  return(results)
  
}


#' @title Modified coefficient of variations (CV)
#' @description Calculate the modified coefficient of variations of a vector, suitable for small sample sizes
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#'
cv = function(x){
  return( (1+1/(4*length(x))) * sd(x)/mean(x) )
}
#
cv.vec = function(x, grp.ids)  ave(x, grp.ids, FUN = cv)
#


#' @title Geometric coefficient of variations (GCV)
#' @description Calculate the geometric coefficient of variations of a vector
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
gcv = function(x){
  return( sqrt( exp(var(log(x))) - 1) )
}
#
gcv.vec = function(x, grp.ids)  ave(x, grp.ids, FUN = gcv)
#
GCV = function(counts, grp.ids=NULL){
  
  # geometric coefficient of variations (GCV)
  # In many applications, it can be assumed that data X are log-normally distributed 
  # (evidenced by the presence of skewness in the sampled data). In such cases, 
  # a more accurate estimate, derived from the properties of the log-normal distribution, is defined as: 
  # Let Y=ln(X) be normally distributed with SD = theta, then c.hat_{v ln} = \sqrt{exp{theta^2}-1}
  
  # this function takes as input a count matrix, and returns GCV for each treatment group (grp) specified for each gene
  # grp = as.factor(c(1,1,1,2,2,2)) indicates 3 reps in each of the two treatments
  
  # problem: produce NaN when a single read count is 0
  
  m = dim(counts)[1]
  n.grp = length(unique(grp.ids))
  gcv.mat = matrix(0, nr=m, nc=n.grp)
  gcv.tmp = t(apply(counts, 1, gcv.vec, grp.ids))
  for (i in 1:m){
    gcv.mat[i, ] = unique(gcv.tmp[i, ])
  }
  return(gcv.mat)
}



#' @title User-defined link function for gamma regression
#' 
#' @description User-defined link function for gamma regression. This is made available to be passed
#' to the \code{family=} argument in the \code{glm} function
#' 
#' @usage
#' family = Gamma(link = vlog())
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
#' @examples
#' 
#' library(NBGOF)
#' set.seed(123)
#' n = 1000                     
#' x = runif(n)
#' sh = 2   
#' vv = vlog()                     
#' y = rgamma(n, scale=vv$linkinv(2+3*x)/sh, shape=sh)
#' fit1 = glm(y~x, family=Gamma(link=vv))                       
#' fit1
#' sum.fit1 = summary(fit1)
#' 1/sum.fit1$dispersion   # compare with sh: close to each other as n increases
#' 
vlog <- function() {
  ## link
  linkfun <- function(y) log(exp(y)-1)
  ## inverse link
  linkinv <- function(eta)  log(exp(eta)+1)
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 1/(exp(-eta) + 1) }
  valideta <- function(eta) TRUE
  link <- "log(exp(y)-1)"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}


