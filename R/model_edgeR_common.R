
#' @title Modeling NB2 Common Dispersion with the Adjusted Profile Likelihood 
#' Estimator (APLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NB regression model with
#' a common dispersion using the adjusted profile likelihood estimator. See details
#' below. The output of this function will be passed to the main GOF function 
#' \code{\link{nb_gof_m}}.
#' 
#' @details This function calls the \code{\link{estimateGLMCommonDisp}} function in
#' \code{edgeR}, with the default estimation method using Cox-Reid adjusted profile likelihood
#' (\code{\link{dispCoxReid}}). See the 
#' \code{\link{estimateGLMCommonDisp}} and \code{\link{glmFit}} functions in the 
#' \code{\link{edgeR}} package for more information.
#' 
#' @usage
#' model_edgeR_common(counts, x, lib.sizes=colSums(counts), design=design)
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
#' @seealso \code{\link{model_edgeR_tagcom}} and \code{\link{model_edgeR_genewise}}
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model_edgeR_common = function(counts, x, lib.sizes=colSums(counts), 
                              design = design){
  
  ## edgeR common dispersion:
  
  stopifnot(design %in% c("single", "multiple", "complex"))
  
  # convert model matrix into group index
  # NOTE: THIS WILL NOT WORK PROPERLY FOR MULTI-GROUP SITUATIONS! (TO DO IN NEXT RELEASE)
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  
  ## the simplest model for a single group (single-intercept model)
  if (design == "single"){    
    design = matrix(as.numeric(as.character(grp.ids)))
    e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
    com.fit = glmFit(d, design, dispersion=e.com$common.dispersion)
    
    # extract quantities:
    mu.hat.m = com.fit$fitted.values  # mu may be close to 0
    phi.hat.m = com.fit$dispersion    # there may be NA's
    v = mu.hat.m + phi.hat.m * mu.hat.m^2
    res.m = (counts - mu.hat.m) / sqrt(v)   # res.m may contain an entire row of NaN
    res.m[is.nan(res.m)] = 0    # replace NaN with 0
    # for 0 reads, sqrt(v) != 0 makes res.m == 0
    # make sure 0/0 = NaN never happens!
    
    # sort res.m with care!
    res.om = t(apply(res.m, 1, sort))  # if NaN is produced in res.m: error when sorting!
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
  
  ## multiple-group case (e.g. two-group comparisons) and complex design
  ## passing x to the design matrix directly
  ## CAUTION: results may not make any sense for complex case!
  if (design == "multiple" | design == "complex"){
    design = x
    e.com = estimateGLMCommonDisp(counts, design, verbose=FALSE)
    com.fit = glmFit(counts, design, dispersion = e.com)
    
    # extract quantities:
    mu.hat.m = com.fit$fitted.values 
    phi.hat.m = com.fit$dispersion  
    v = mu.hat.m + phi.hat.m * mu.hat.m^2
    res.m = (counts - mu.hat.m) / sqrt(v)
    if (typeof(res.m) == "list"){
      dimension = dim(res.m)
      unlist.res.m = unlist(res.m)
      unlist.res.m[ is.nan(unlist.res.m) ] = 0
      res.m = matrix( unlist.res.m, dimension )
    }
    else 
      res.m[is.nan(res.m)] = 0
    
    # sort res.m with care!
    # res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
    res.om = t(apply(res.m, 1, sort))  
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
  
}