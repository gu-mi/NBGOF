
# dispersion estimation: trended dispersion in edgeR

################################################################################
#' @title Modeling NB trended dispersion model with the adjusted profile likelihood
#' estimator (APLE) on original and simulated datasets
#' 
#' @description This function is designed to fit an NB regression model with
#' trended dispersions using the adjusted profile likelihood estimator. In edgeR, 
#' this function assumes the dispersion \eqn{\phi_i} satisfies 
#' \eqn{\phi_i=s(\bar{\mu}_{i\cdot})}, where \eqn{s(\cdot)} is a smooth function of 
#' each gene's average read counts across samples. A variety of non-parametric 
#' approaches can be used by fitting loess or spline curves on binned genes, or 
#' using locally weighted APL. See details  below. The output of this function will 
#' be passed to the main GOF function.
#' 
#' @details Under the NB model ... See the 
#' \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTagwiseDisp}} and
#' \code{\link{glmFit}} functions in the \code{\link{edgeR}} package
#' for more information.
#' 
#' @usage
#' model_edgeR_trended(counts, x, lib.sizes=colSums(counts))
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @seealso \code{\link{model_edgeR_genewise}} for the non-shrinkage version
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_edgeR_trended <- function(counts, x, lib.sizes=colSums(counts)){

  ## edgeR tagwise dispersion:
  
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  design = model.matrix(~grp.ids, data=d$samples)
  e.com = estimateGLMCommonDisp(d, design, verbose=TRUE)
  e.trd = estimateGLMTrendedDisp(e.com, design)
  trd.fit = glmFit(d, design, dispersion=e.trd$trended.dispersion)
  
  # extract quantities:
  mu.hat.m = trd.fit$fitted.values   # mu may be close to 0
  phi.hat.m = trd.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
  res.m = (counts - mu.hat.m) / sqrt(v) 
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
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