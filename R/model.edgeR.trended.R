
#' @title Modeling NB Trended (Non-parametric) Dispersion with the Adjusted Profile
#' Likelihood estimator (APLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NB regression model with
#' trended (non-parametric) dispersions using the adjusted profile likelihood estimator. 
#' In \code{edgeR}, this function assumes the dispersion \eqn{\phi_i} satisfies 
#' \eqn{\phi_i=s(\bar{\mu}_{i\cdot})}, where \eqn{s(\cdot)} is a smooth function of 
#' each gene's average read counts across samples. A variety of non-parametric 
#' approaches can be used by fitting loess or spline curves on binned genes, or 
#' using locally weighted APL. See details below. The output of this function will 
#' be passed to the main GOF function \code{\link{nb.gof.m}}.
#' 
#' @details In this trended non-parametric model, \eqn{\phi_{ij}} is estimated in a 
#' first step as a smooth function of 
#' \eqn{\log(\hat{\phi}_{ij})} on \eqn{\log(\hat{\mu}_{ij})}, and then treated as known 
#' in the second step of regression coefficient inference. See the 
#' \code{\link{estimateGLMTrendedDisp}} and
#' \code{\link{glmFit}} functions in the \code{\link{edgeR}} package
#' for more information.
#' 
#' @usage
#' model.edgeR.trended(counts, x, lib.sizes=colSums(counts), min.n=min.n, method=method)
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
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model.edgeR.trended = function(counts, x, lib.sizes=colSums(counts), min.n=min.n, method=method){
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  
  ## edgeR trended (non-parametric) dispersion:

  method = ifelse(test = is.null(method), "auto", method)
  stopifnot(method %in% c("auto", "bin.spline", "bin.loess", "power", "spline"))
  
  y.dge = DGEList(counts=counts)
  y.dge$offset = log(lib.sizes)  
  e.trd = estimateGLMTrendedDisp(y.dge, design=x, min.n=min.n, method=method)
  trd.fit = glmFit(y=counts, design=x, dispersion=e.trd$trended.dispersion)
  
  # extract quantities:
  mu.hat.m = trd.fit$fitted.values   # mu may be close to 0
  phi.hat.m = trd.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2
  res.m = as.matrix((counts - mu.hat.m) / sqrt(v))
  
  # make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
  res.m[ is.nan(res.m) ] = 0
  res.m[ is.infinite(res.m) ] = 0
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_trd_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
  )
  return(model_trd_m_obj)
}



