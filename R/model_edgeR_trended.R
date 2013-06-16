
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
#' be passed to the main GOF function \code{\link{nb_gof_m}}.
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
#' model_edgeR_trended(counts, x, lib.sizes=colSums(counts), min.n=min.n, design=design)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param min.n minimim number of genes in a bin. Default is 100. See \code{\link{dispBinTrend}}
#' for details (lower-level function of \code{\link{estimateGLMTrendedDisp}}).
#' @param design specifications of testing (1) a single-group model (\code{single}); (2)
#' a multiple-group model (\code{multiple}); (3) complex design with interactions
#' (\code{complex}). \strong{The codes are only tested for single- and multiple-group cases, and
#' the modeling of dispersions in the multiple-group case is still under consideration!
#' For the complex design with interactions, though we can pass the design matrix directly
#' in edgeR, the results may not make any sense!}
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model_edgeR_trended = function(counts, x, lib.sizes=colSums(counts), min.n=min.n, design=design){
  
  ## edgeR trended (non-parametric) dispersion:
  
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  
  ## the simplest model for a single group (single-intercept model)
  if (design == "simple"){
      if (length(unique(grp.ids)) == 1){   
      design = matrix(as.numeric(as.character(grp.ids)))
      e.trd = estimateGLMTrendedDisp(d, design, min.n=min.n)
      trd.fit = glmFit(d, design, dispersion=e.trd$trended.dispersion)
      
      # extract quantities:
      mu.hat.m = trd.fit$fitted.values   # mu may be close to 0
      phi.hat.m = trd.fit$dispersion     # there may be NA's
      v = mu.hat.m + phi.hat.m * mu.hat.m^2
      res.m = (counts - mu.hat.m) / sqrt(v)
      res.m[is.nan(res.m)] = 0
      
      # sort res.m with care!
      res.om = t(apply(res.m, 1, sort))
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
  }
    
  ## multiple-group case (e.g. two-group comparisons)
  if (design == "multiple") { 
      design = model.matrix(~grp.ids, data=d$samples)
      e.trd = estimateGLMTrendedDisp(d, design, min.n=min.n)
      trd.fit = glmFit(d, design, dispersion=e.trd$trended.dispersion)
      
      # extract quantities:
      mu.hat.m = trd.fit$fitted.values   # mu may be close to 0
      phi.hat.m = trd.fit$dispersion     # there may be NA's
      v = mu.hat.m + phi.hat.m * mu.hat.m^2
      res.m = (counts - mu.hat.m) / sqrt(v)
      res.m[is.nan(res.m)] = 0
      
      # sort res.m with care!
      res.om = t(apply(res.m, 1, sort))
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
  
  ## complex design: passing x to the design matrix directly
  ## CAUTION: results may not make any sense!
  if (design == "complex"){
    design = x
    e.trd = estimateGLMTrendedDisp(d, design, min.n=min.n)
    trd.fit = glmFit(d, design, dispersion=e.trd$trended.dispersion)
    
    # extract quantities:
    mu.hat.m = trd.fit$fitted.values   # mu may be close to 0
    phi.hat.m = trd.fit$dispersion     # there may be NA's
    v = mu.hat.m + phi.hat.m * mu.hat.m^2
    res.m = (counts - mu.hat.m) / sqrt(v)
    res.m[is.nan(res.m)] = 0
    
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
}

