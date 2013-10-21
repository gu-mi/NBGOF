
#' @title Modeling NB Tagwise-Common Dispersion with the Adjusted Profile Likelihood
#' Estimator (APLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NB regression model with
#' tagwise dispersions using the adjusted profile likelihood estimator (APLE). In \code{edgeR}, 
#' this tagwise dispersion model lies between two extreme scenarios: entirely 
#' individual \eqn{\phi_i} for each gene and entirely shared values (i.e. common dispersion).
#' The function \code{estimateGLMTagwiseDisp} in \code{edgeR} shrinks the dispersion towards 
#' a common dispersion or a global dispersion trend. In our implementation, the function
#' \code{model.edgeR.tagcom} 
#' shrinks the tagwise dispersions towards a common dispersion. See details below. The output 
#' of this function will be passed to the main GOF function \code{\link{nb.gof.m}}.
#' 
#' @details In this tagwise-common dispersion model, the dispersion parameter \eqn{\phi_i} is
#' estimated by maximizing a penalized log-likelihood \eqn{APL_g(\phi_g)} plus a weighted shared 
#' log-likelihood. The weight denoted by \eqn{G_0} controls the level of shrinkage applied to 
#' purely genewise dispersions towards a common dispersion. In the paper of 
#' McCarthy et al. (2012),
#' \eqn{G_0=20/df} performs well over real RNA-Seq datasets. Here the numerator 20 is the prior
#' degrees of freedom, which can be specified by the argument \code{prior.df} in the lower-level
#' function \code{\link{dispCoxReidInterpolateTagwise}} in \code{edgeR}. The denominator "df" is
#' the residual degrees of freedom (number of libraries minus the number of distinct treatment
#' groups). We leave the default value of \code{prior.df} at 10. See the 
#' \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTagwiseDisp}}, 
#' \code{\link{dispCoxReidInterpolateTagwise}} and \code{\link{glmFit}} functions in the 
#' \code{edgeR} package for more details.
#' 
#' @usage
#' model.edgeR.tagcom(counts, x, lib.sizes=colSums(counts), prior.df = prior.df, method=method)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param prior.df prior degrees of freedom to control the level of shrinkage (default is 10). See 
#' \code{\link{estimateGLMTagwiseDisp}} for more details.
#' This argument is controlled by a higher level function \code{\link{nb.gof.m}}.
#' @param method method for estimating the common dispersion, including "CoxReid", "Pearson" and "deviance". If NULL, then
#' use the "CoxReid" method. See \code{\link{estimateGLMCommonDisp}} for more details.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.m}} function.
#' 
#' @seealso \code{\link{model.edgeR.genewise}} for the non-shrinkage, entirely genewise model,
#' and \code{\link{model.edgeR.tagtrd}} for the tagwise model that shrinks towards the trended
#' dispersion.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model.edgeR.tagcom = function(counts, x, lib.sizes=colSums(counts), prior.df = prior.df, method=method){
  
  ## edgeR tagwise-common dispersion:
  
  method = ifelse(test = is.null(method), "CoxReid", method)
  stopifnot(method %in% c("CoxReid", "Pearson", "deviance"))
  
  e.com = estimateGLMCommonDisp(y=counts, design=x, verbose=FALSE, method=method)
  e.tgc = estimateGLMTagwiseDisp(y=counts, design=x, offset=log(lib.sizes), dispersion=e.com, prior.df = prior.df, trend=FALSE)
  # e.tag: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
  tgc.fit = glmFit(y=counts, design=x, dispersion=e.tgc$tagwise.dispersion)
  
  # extract quantities:
  mu.hat.m = tgc.fit$fitted.values   # mu may be close to 0
  phi.hat.m = tgc.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2
  res.m = (counts - mu.hat.m) / sqrt(v)
  
  # make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
  res.m[ is.nan(res.m) ] = 0
  res.m[ is.infinite(res.m) ] = 0
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort))
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_tgc_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
  )
  return(model_tgc_m_obj)
}