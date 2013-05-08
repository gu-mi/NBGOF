
#' @title Modeling NB2 common dispersion model with the adjusted profile likelihood 
#' estimator (APLE) on original and simulated datasets
#' 
#' @description This function is designed to fit an NB regression model with
#' a common dispersion using the adjusted profile likelihood estimator. See details
#' below. The output of this function will be passed to the main GOF function 
#' \code{\link{nb_gof_m}}.
#' 
#' @details Under the NB model ... See the 
#' \code{\link{estimateGLMCommonDisp}} and \code{\link{glmFit}} functions in the 
#' \code{\link{edgeR}} package for more information.
#' 
#' @usage
#' model_edgeR_common(counts, x, lib.sizes=colSums(counts))
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @seealso \code{\link{model_edgeR_tagwise}} and \code{\link{model_edgeR_trended}}
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_edgeR_common <- function(counts, x, lib.sizes=colSums(counts), design = design){
  
  ## edgeR common dispersion:
  
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
      com.fit = glmFit(d, design, dispersion=e.com$common.dispersion)
      
      # extract quantities:
      mu.hat.m = com.fit$fitted.values  # mu may be close to 0
      phi.hat.m = com.fit$dispersion    # there may be NA's
      v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14  # add epsilon to avoid v == 0
      res.m = (counts - mu.hat.m) / sqrt(v)   
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
    
    # multiple-group case
    else {    
      design = model.matrix(~grp.ids, data=d$samples)
      e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
      com.fit = glmFit(d, design, dispersion=e.com$common.dispersion)
      
      # extract quantities:
      mu.hat.m = com.fit$fitted.values 
      phi.hat.m = com.fit$dispersion  
      v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
      res.m = (counts - mu.hat.m) / sqrt(v)    
      
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
  }
  
  else if (design == "complex"){
    design = x
    e.com = estimateGLMCommonDisp(d, design, verbose=FALSE)
    com.fit = glmFit(d, design, dispersion=e.com$common.dispersion)
    
    # extract quantities:
    mu.hat.m = com.fit$fitted.values 
    phi.hat.m = com.fit$dispersion  
    v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
    res.m = (counts - mu.hat.m) / sqrt(v)    
    
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
  
}