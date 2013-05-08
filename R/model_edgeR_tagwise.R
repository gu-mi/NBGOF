
#' @title Modeling NB tagwise dispersion model with the adjusted profile likelihood
#' estimator (APLE) on original and simulated datasets
#' 
#' @description This function is designed to fit an NB regression model with
#' genewise dispersions using the adjusted profile likelihood estimator. In \code{edgeR}, 
#' this tagwise dispersion model lies between two extreme scenarios: entirely 
#' individual \eqn{\phi_i} for each gene and entirely shared values (i.e. common dispersion).
#' This function shrinks the dispersion towards a common dispersion or a global dispersion 
#' trend. See details below. The output of this function will be passed to the main GOF function
#' \code{\link{nb_gof_m}}.
#' 
#' @details In this tagwise dispersion model, the dispersion parameter $\phi_i$ is estimated by
#' maximizing a penalized log-likelihood $APL_g(\phi_g)$ plus a weighted shared log-likelihood.
#' The weight denoted by $G_0$ controls the level of shrinkage of purely genewise dispersions
#' towards a common dispersion. In the paper of McCarthy $et al.$ (2012), $G_0=20/df$ performs
#' well over real RNA-Seq datasets. Here the numerator 20 is the prior degrees of freedom that can
#' be specified (\code{prior.df}) in the lower-level function 
#' \code{\link{dispCoxReidInterpolateTagwise}} in \code{\link{edgeR}}. The denominator "df" is the
#' residual degrees of freedom (number of libraries minus the number of distinct treatment groups).
#' We leave the default value of \code{prior.df} at 10. See the 
#' \code{\link{estimateGLMCommonDisp}}, \code{\link{estimateGLMTagwiseDisp}}, 
#' \code{\link{dispCoxReidInterpolateTagwise}} and \code{\link{glmFit}} functions in the 
#' \code{\link{edgeR}} package for more information.
#' 
#' @usage
#' model_edgeR_tagwise(counts, x, lib.sizes=colSums(counts), prior.df = prior.df)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param prior.df prior degrees of freedom to control the level of shrinkage (default is 10)
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @seealso \code{\link{model_edgeR_genewise}} for the non-shrinkage, entirely genewise version
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_edgeR_tagwise <- function(counts, x, lib.sizes=colSums(counts), prior.df = prior.df,
                                design = design){
  
  ## edgeR genewise dispersion:
  
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  
  if (design == "simple"){
    
    # include fitting a single-intercept model, so separate into two parts:
    
    # single-group case
    if (length(unique(grp.ids)) == 1){   
      design = matrix(as.numeric(as.character(grp.ids)))
      e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
      e.tag = estimateGLMTagwiseDisp(e.com, design, prior.df = prior.df)
      # e.tag: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
      tag.fit = glmFit(d, design, dispersion=e.tag$tagwise.dispersion)
      
      # extract quantities:
      mu.hat.m = tag.fit$fitted.values   # mu may be close to 0
      phi.hat.m = tag.fit$dispersion     # there may be NA's
      v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
      res.m = (counts - mu.hat.m) / sqrt(v) 
      
      # sort res.m with care!
      res.om = t(apply(res.m, 1, sort))
      ord.res.v = as.vector(t(res.om))
      
      # save as a list
      model_tag_m_obj = list(mu.hat.mat = mu.hat.m,
                             res.mat = res.m,
                             res.omat = res.om,
                             ord.res.vec = ord.res.v,
                             phi.hat.mat = phi.hat.m
      )
      return(model_tag_m_obj)
    }
    
    # multiple-group case
    else { 
      design = model.matrix(~grp.ids, data=d$samples)
      e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
      e.tag = estimateGLMTagwiseDisp(e.com, design, prior.df = prior.df)
      # e.tag: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
      tag.fit = glmFit(d, design, dispersion=e.tag$tagwise.dispersion)
      
      # extract quantities:
      mu.hat.m = tag.fit$fitted.values   # mu may be close to 0
      phi.hat.m = tag.fit$dispersion     # there may be NA's
      v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
      res.m = (counts - mu.hat.m) / sqrt(v) 
      
      # sort res.m with care!
      res.om = t(apply(res.m, 1, sort.vec, grp.ids))
      ord.res.v = as.vector(t(res.om))
      
      # save as a list
      model_tag_m_obj = list(mu.hat.mat = mu.hat.m,
                             res.mat = res.m,
                             res.omat = res.om,
                             ord.res.vec = ord.res.v,
                             phi.hat.mat = phi.hat.m
      )
      return(model_tag_m_obj)
    } 
  }
  
  else if (design == "complex"){
    design = x
    e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
    e.tag = estimateGLMTagwiseDisp(e.com, design, prior.df = prior.df)
    # e.tag: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
    tag.fit = glmFit(d, design, dispersion=e.tag$tagwise.dispersion)
    
    # extract quantities:
    mu.hat.m = tag.fit$fitted.values   # mu may be close to 0
    phi.hat.m = tag.fit$dispersion     # there may be NA's
    v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
    res.m = (counts - mu.hat.m) / sqrt(v) 
    
    # sort res.m with care!
    res.om = t(apply(res.m, 1, sort.vec, grp.ids))
    ord.res.v = as.vector(t(res.om))
    
    # save as a list
    model_tag_m_obj = list(mu.hat.mat = mu.hat.m,
                           res.mat = res.m,
                           res.omat = res.om,
                           ord.res.vec = ord.res.v,
                           phi.hat.mat = phi.hat.m
    )
    return(model_tag_m_obj)
  }

}