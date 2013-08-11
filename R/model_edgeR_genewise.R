
#' @title Modeling NB Genewise Dispersion with the Adjusted Profile Likelihood
#' Estimator (APLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NB regression model with
#' genewise dispersions using the adjusted profile likelihood estimator. This function
#' differs from the \code{\link{model_edgeR_tagcom}} and \code{\link{model_edgeR_tagtrd}}
#' functions: the tagwise model shrinks the dispersion towards a common dispersion or a global
#' dispersion trend, while the genewise model (implemented by this function) applies 
#' \strong{no shrinkage}. This model is \emph{not} recommended, and is used for diagnostics
#' only. See details below. The output of this function 
#' will be passed to the main GOF function \code{\link{nb_gof_m}}.
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
#' model_edgeR_genewise(counts, x, lib.sizes=colSums(counts), design=design)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param design specifications of testing (1) a single-group model (\code{single}); (2)
#' a multiple-group model (\code{multiple}); (3) complex design with interactions
#' (\code{complex}). \strong{The codes are only tested for single- and multiple-group cases, and
#' the modeling of dispersions in the multiple-group case is still under consideration!
#' For the complex design with interactions, though we can pass the design matrix directly
#' in edgeR, the results may not make any sense!}
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @seealso \code{\link{model_edgeR_tagcom}} and \code{\link{model_edgeR_tagtrd}} for the
#' shrinkage versions.
#'  
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model_edgeR_genewise = function(counts, x, lib.sizes=colSums(counts), design=design){
  
  ## edgeR genewise dispersion:
  
  stopifnot(design %in% c("single", "multiple", "complex"))
  
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  
  ## the simplest model for a single group (single-intercept model)
  if (design == "single"){
    design = matrix(as.numeric(as.character(grp.ids)))
    e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
    e.gen = estimateGLMTagwiseDisp(e.com, design, prior.df=0)  # prior.df = 0
    # e.gen: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
    gen.fit = glmFit(d, design, dispersion=e.gen$tagwise.dispersion)
    
    # extract quantities:
    mu.hat.m = gen.fit$fitted.values   # mu may be close to 0
    phi.hat.m = gen.fit$dispersion     # there may be NA's
    v = mu.hat.m + phi.hat.m * mu.hat.m^2
    res.m = (counts - mu.hat.m) / sqrt(v)
    res.m[is.nan(res.m)] = 0
    
    # sort res.m with care!
    res.om = t(apply(res.m, 1, sort))
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
    
  ## multiple-group case (e.g. two-group comparisons)
  if (design == "multiple-group") { 
      design = model.matrix(~grp.ids, data=d$samples)
      e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
      e.gen = estimateGLMTagwiseDisp(e.com, design, prior.df=0)  # prior.df = 0
      # e.gen: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
      gen.fit = glmFit(d, design, dispersion=e.gen$tagwise.dispersion)
      
      # extract quantities:
      mu.hat.m = gen.fit$fitted.values   # mu may be close to 0
      phi.hat.m = gen.fit$dispersion     # there may be NA's
      v = mu.hat.m + phi.hat.m * mu.hat.m^2
      res.m = (counts - mu.hat.m) / sqrt(v)
      res.m[is.nan(res.m)] = 0
      
      # sort res.m with care!
      res.om = t(apply(res.m, 1, sort))
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
  
  ## complex design: passing x to the design matrix directly
  ## CAUTION: results may not make any sense!
  if (design == "complex"){
    design = x
    e.com = estimateGLMCommonDisp(counts, design, verbose=FALSE)
    e.gen = estimateGLMTagwiseDisp(e.com, design, prior.df=0)  # prior.df = 0
    gen.fit = glmFit(counts, design, dispersion=e.gen$tagwise.dispersion)
    
    # extract quantities:
    mu.hat.m = gen.fit$fitted.values   # mu may be close to 0
    phi.hat.m = gen.fit$dispersion     # there may be NA's
    v = mu.hat.m + phi.hat.m * mu.hat.m^2
    res.m = (counts - mu.hat.m) / sqrt(v)
    res.m[is.nan(res.m)] = 0
    
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
}

  
