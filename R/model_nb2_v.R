
## For NB2 model fitting on the original & simulated datasets

# the same function, nb.regression.1(), is used in NB2 fitting for matrix response
# cf. model_nb2_m.R

################################################################################
#' @title Modeling NB2 regression model with MLE on original and simulated 
#' datasets
#' 
#' @description This function is designed to fit an NB2 regression model. The output of
#' this function will be passed to the main GOF function.
#' 
#' @details Details here
#' 
#' @usage
#' model_nb2_v(y, x, lib.sizes=NULL)
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
model_nb2_v <- function(y, x, lib.sizes=NULL){
  
  n = length(y)
  p = dim(x)[2]
  lib.sizes = ifelse(rep(is.null(lib.sizes),n), rep(1,n), lib.sizes)
  
  # preconditions
  stopifnot(n == dim(x)[1], n == length(lib.sizes))
  
  # fit NB2 model using MLE
  nb2.fit = nb.regression.1(y=y, s=lib.sizes, x=x, beta=NA)
  mu.hat.v = nb2.fit$mu
  phi = nb2.fit$phi
  v = nb2.fit$v
  res.v = (y - mu.hat.v)/sqrt(v)
  
  # save as a list
  model_nb2_v_obj = list(mu.hat.v = mu.hat.v,
                         res.vec = res.v,
                         phi = phi
                         )
  return(model_nb2_v_obj)
}
