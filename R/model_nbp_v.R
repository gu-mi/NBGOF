
## For NBP model fitting on the original & simulated dataset
## Now include both MLE and APLE estimation methods

################################################################################
#' @title Modeling NBP regression model with MLE or APLE on original and simulated 
#' datasets
#' 
#' @description This function is designed to fit an NBP regression model. The output of
#' this function will be passed to the main GOF function.
#' 
#' @details Details here
#' 
#' @usage
#' model_nbp_v(y, x, lib.sizes=NULL, est.method="MLE")
#' 
#' @param y an n-by-1 vector of non-negative integers. For a typical RNA-Seq experiment, 
#' this may represent the read counts for a single gene 
#' @param x an n-by-p design matrix. If an intercept is desired in the model, you need to specify
#' the first column of \code{x} as a vector of 1.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is 1 for all samples
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_v}} function.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_nbp_v <- function(y, x, lib.sizes=NULL, est.method = "MLE"){
  
  n = length(y)
  p = dim(x)[2]
  lib.sizes = ifelse(rep(is.null(lib.sizes),n), rep(1,n), lib.sizes)
  
  # preconditions
  stopifnot(n == dim(x)[1], n == length(lib.sizes), est.method %in% c("MLE", "APLE"))
  
  if (est.method == "MLE"){
    # fit NBP model using MLE: glm.nbp.1.MLE()
    fit = glm.nbp.1.MLE(y=y, s=lib.sizes, x=x, print.level=0)
    mu = fit$mu
    phi = fit$phi
    res.v = fit$p.res
  }

  else if (est.method == "APLE"){
    # fit NBP model using APLE: glm.nbp.1()
    fit = glm.nbp.1(y=y, s=lib.sizes, x=x, print.level=0)
    mu = fit$mu
    phi = fit$phi
    res.v = fit$p.res
  }

  # save as a list
  model_nbp_v_obj = list(mu.hat.v = mu,
                         res.vec = res.v,
                         phi = phi
                         )
  return(model_nbp_v_obj)
}
