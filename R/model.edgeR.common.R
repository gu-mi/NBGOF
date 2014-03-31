
#' @title Modeling NB2 Common Dispersion with the Adjusted Profile Likelihood 
#' Estimator (APLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NB regression model with
#' a common dispersion using the adjusted profile likelihood estimator. See details
#' below. The output of this function will be passed to the main GOF function 
#' \code{\link{nb.gof.m}}.
#' 
#' @details This function calls the \code{\link{estimateGLMCommonDisp}} function in
#' \code{edgeR}, with the default estimation method using Cox-Reid adjusted profile likelihood
#' (\code{\link{dispCoxReid}}). See the 
#' \code{\link{estimateGLMCommonDisp}} and \code{\link{glmFit}} functions in the 
#' \code{\link{edgeR}} package for more information.
#' 
#' @usage
#' model.edgeR.common(counts, x, lib.sizes=colSums(counts), method=method)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param method method for estimating the common dispersion, including "CoxReid", "Pearson" and "deviance". If NULL, then
#' use the "CoxReid" method. See \code{\link{estimateGLMCommonDisp}} for more details.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.m}} function.
#' 
#' @seealso \code{\link{model.edgeR.tagcom}} and \code{\link{model.edgeR.genewise}}
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model.edgeR.common = function(counts, x, lib.sizes=colSums(counts), method=method){
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  
  ## edgeR common dispersion:
  
  method = ifelse(test = is.null(method), "CoxReid", method)
  stopifnot(method %in% c("CoxReid", "Pearson", "deviance"))
  
  e.com = estimateGLMCommonDisp(y=counts, design=x, method=method, verbose=FALSE)
  com.fit = glmFit(y=counts, design=x, dispersion=e.com)
  
  # extract quantities:
  mu.hat.m = com.fit$fitted.values  # mu may be close to 0
  phi.hat.m = com.fit$dispersion    # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2
  res.m = as.matrix((counts - mu.hat.m) / sqrt(v))   # res.m may contain an entire row of NaN
  
  # make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
  # make sure res.m is a matrix, not a list
  res.m[ is.nan(res.m) ] = 0
  res.m[ is.infinite(res.m) ] = 0
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_com_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
  )
  return(model_com_m_obj) 
}

