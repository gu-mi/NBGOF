
## For NB2 model fitting on the original & simulated datasets
# dispersion estimation: tagwise dispersion in edgeR (default with shrinkage)

# update: 2013/01/11

################################################################################
#' @title Modeling NB tagwise dispersion model with the adjusted profile likelihood
#' estimator (APLE) on original and simulated datasets
#' 
#' @description This function is designed to fit an NB regression model with
#' genewise dispersions using the adjusted profile likelihood estimator. In edgeR, 
#' this tagwise dispersion model actually lie between two extreme scenarios: entirely 
#' individual \eqn{\phi_i} for each gene and entirely shared values. This function 
#' shrinkages the dispersion towards the global dispersion trend. See details
#' below. The output of this function will be passed to the main GOF function.
#' 
#' @details Under the NB model ... See the 
#' \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTagwiseDisp}} and
#' \code{\link{glmFit}} functions in the \code{\link{edgeR}} package
#' for more information.
#' 
#' @usage
#' model_edgeR_tagwise(counts, x, lib.sizes=colSums(counts))
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
model_edgeR_tagwise <- function(counts, x, lib.sizes=colSums(counts)){
  
#   nr = dim(counts)[1]
#   nc = dim(counts)[2]
#   n = nr * nc
  
  # initializations
#   mu.hat.m = matrix(0, nr = nr, nc = nc)
#   phi.hat.m = matrix(0, nr = nr, nc = nc)
#   res.m = matrix(0, nr = nr, nc = nc)
  
  ## edgeR tagwise dispersion:
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  design = model.matrix(~grp.ids, data=d$samples)
  e.com = estimateGLMCommonDisp(d, design, verbose=TRUE)
  e.tag = estimateGLMTagwiseDisp(e.com, design)
  # e.tag: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
  tag.fit = glmFit(d, design, dispersion=e.tag$tagwise.dispersion)
  
  # extract quantities:
  mu.hat.m = tag.fit$fitted.values   # mu may be close to 0
  phi.hat.m = tag.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
  res.m = (counts - mu.hat.m) / sqrt(v) 
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_tag_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
                         #fit.tag = tag.fit
  )
  return(model_tag_m_obj)
}