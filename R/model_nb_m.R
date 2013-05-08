
#' @title Modeling NB genewise dispersion model with MLE on original and simulated 
#' datasets
#' 
#' @description This function is designed to fit an NB regression model with
#' genewise dispersions using the maximum likelihood estimator. The output of
#' this function will be passed to the main GOF function.
#' 
#' @details Details here
#' 
#' @usage
#' model_nb_m(counts, x, lib.sizes=colSums(counts))
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_nb_m <- function(counts, x, lib.sizes=colSums(counts)){
  
  nc = dim(counts)[2]  
  
  # preconditions
  stopifnot(is.matrix(x), nc == dim(x)[1])
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  
  # include fitting a single-intercept model, so separate into two parts:
  
  # single-group case
  if (length(unique(grp.ids)) == 1){   
    # NB2 fit and extract quantities
    nb.fit = nb.regression.1(y=counts, s=lib.sizes, x=x, beta=NA)
    mu.hat.m = nb.fit$mu
    v.hat.m = nb.fit$v
    phi.hat.m = nb.fit$phi
    res.m = (counts - mu.hat.m) / sqrt(v.hat.m) 
    
    # sort res.m with care!
    res.om = t(apply(res.m, 1, sort))
    ord.res.v = as.vector(t(res.om))   # the "V" vector of block-wise-ordered residuals
    
    # save as a list
    model_nb_m_obj = list(mu.hat.mat = mu.hat.m,
                          res.mat = res.m,
                          res.omat = res.om,
                          ord.res.vec = ord.res.v,
                          phi.hat.mat = phi.hat.m
    )
    return(model_nb_m_obj)
  }
  
  # multiple-group case
  else { 
    # NB2 fit and extract quantities
    nb.fit = nb.regression.1(y=counts, s=lib.sizes, x=x, beta=NA)
    mu.hat.m = nb.fit$mu
    v.hat.m = nb.fit$v
    phi.hat.m = nb.fit$phi
    res.m = (counts - mu.hat.m) / sqrt(v.hat.m) 
    
    # sort res.m with care!
    res.om = t(apply(res.m, 1, sort.vec, grp.ids))
    ord.res.v = as.vector(t(res.om))   # the "V" vector of block-wise-ordered residuals
    
    # save as a list
    model_nb_m_obj = list(mu.hat.mat = mu.hat.m,
                          res.mat = res.m,
                          res.omat = res.om,
                          ord.res.vec = ord.res.v,
                          phi.hat.mat = phi.hat.m
    )
    return(model_nb_m_obj)
  }
    
}
