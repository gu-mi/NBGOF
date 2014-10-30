
#' @title Modeling NB Genewise Dispersion with the Adjusted Profile Likelihood
#' Estimator (APLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NB regression model with
#' genewise dispersions using the adjusted profile likelihood estimator. This function
#' differs from the \code{\link{model.edgeR.tagcom}} and \code{\link{model.edgeR.tagtrd}}
#' functions: the tagwise model shrinks the dispersion towards a common dispersion or a global
#' dispersion trend, while the genewise model (implemented by this function) applies 
#' \strong{no shrinkage}. This model is \emph{not} recommended, and is used for diagnostics
#' only. See details below. The output of this function 
#' will be passed to the main GOF function \code{\link{nb.gof.m}}.
#' 
#' @details In this genewise dispersion model, the dispersion parameter \eqn{\phi_i} is
#' estimated by maximizing a penalized log-likelihood \eqn{APL_g(\phi_g)} plus a weighted shared 
#' log-likelihood. The weight denoted by \eqn{G_0} controls the level of shrinkage applied to 
#' purely genewise dispersions towards a common dispersion. In the paper of 
#' McCarthy et al. (2012),
#' \eqn{G_0=20/df} performs well over real RNA-Seq datasets. Here the numerator 20 is the prior
#' degrees of freedom, which can be specified by the argument \code{prior.df} in the lower-level
#' function \code{\link{dispCoxReidInterpolateTagwise}} in \code{edgeR}. The denominator "df" is
#' the residual degrees of freedom (number of libraries minus the number of distinct treatment
#' groups). We set the value of \code{prior.df} equals 0, so that no shrinkage is applied. 
#' See the 
#' \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTagwiseDisp}}, 
#' \code{\link{dispCoxReidInterpolateTagwise}} and \code{\link{glmFit}} functions in the 
#' \code{edgeR} package for more details.
#' 
#' @usage
#' model.edgeR.genewise(counts, x, lib.sizes=colSums(counts), min.n=min.n, method=method)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param min.n minimim number of genes in a bin. Default is 100. See \code{\link{dispBinTrend}}
#' for details (lower-level function of \code{\link{estimateGLMTrendedDisp}}).
#' @param method method for estimating the trended dispersion, including "auto", "bin.spline", "bin.loess", "power" and "spline". 
#' If NULL, then the "auto" method. Normally the number of tags analyzed is greater than 200, so the "bin.spline" method is used which
#' calls the \code{\link{dispBinTrend}} function. See \code{\link{estimateGLMTrendedDisp}} for more details.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.m}} function.
#' 
#' @seealso \code{\link{model.edgeR.tagcom}} and \code{\link{model.edgeR.tagtrd}} for the
#' shrinkage versions.
#'  
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
#' @export
#' 
model.edgeR.genewise = function(counts, x, lib.sizes=colSums(counts), min.n = min.n, method=method){
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  
  ## edgeR genewise dispersion:
  
  method = ifelse(test = is.null(method), "auto", method)
  stopifnot(method %in% c("auto", "bin.spline", "bin.loess", "power", "spline"))
  
  y.dge = DGEList(counts = y, group = grp.ids)
  y.dge = calcNormFactors(y.dge, method="RLE")
  y.dge = estimateGLMTrendedDisp(y.dge, design=x, min.n=min.n, verbose=FALSE)
  e.gen = estimateGLMTagwiseDisp(y.dge, design=x, prior.df=0, trend=TRUE)  # prior.df = 0
  # e.gen: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
  gen.fit = glmFit(e.gen, design=x, dispersion=e.gen$tagwise.dispersion)
  
  
#   y.dge = DGEList(counts=counts)
#   y.dge$offset = log(lib.sizes)  
#   y.dge = estimateGLMTrendedDisp(y.dge, design=x, min.n=min.n, method=method)
#   e.gen = estimateGLMTagwiseDisp(y.dge, design=x, dispersion=y.dge$trended.dispersion, prior.df = 0, trend=TRUE)  # prior.df = 0
#   gen.fit = glmFit(y=counts, design=x, dispersion=e.gen$tagwise.dispersion)
  
  # extract quantities:
  mu.hat.m = gen.fit$fitted.values   # mu may be close to 0
  phi.hat.m = gen.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2
  res.m = as.matrix((counts - mu.hat.m) / sqrt(v))
  
  # make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
  res.m[ is.nan(res.m) ] = 0
  res.m[ is.infinite(res.m) ] = 0
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_gen_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
  )
  return(model_gen_m_obj)
}




