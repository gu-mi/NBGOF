
#' @title Modeling Poisson Regression Model on Original and Simulated Datasets
#' 
#' @description This function is designed to fit a Poisson regression model. The output of
#' this function will be passed to the main GOF function.
#' 
#' @details The \code{glm} function with \code{family=poisson} is used for Poisson model fitting.
#' 
#' @usage
#' model.poi.v(y, x, lib.sizes=NULL)
#' 
#' @param y an n-by-1 vector of non-negative integers. For a typical RNA-Seq experiment, 
#' this may represent the read counts for a single gene across n samples.
#' @param x an n-by-p design matrix. For Poisson model fitting, we used the \link{glm} function,
#' so if an intercept is desired, there is no need to include the first column of 1.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is 0 for all samples.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.v}} function.
#' 
#' @author Gu Mi <neo.migu@gmail.com>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model.poi.v <- function(y, x, lib.sizes=NULL){
  
  n = length(y)
  p = dim(x)[2]
  lib.sizes = ifelse(rep(is.null(lib.sizes),n), rep(0,n), lib.sizes)
  
  # preconditions
  stopifnot(n == dim(x)[1], n == length(lib.sizes))
  
  # fit Poisson model
  poi.fit = glm(y ~ x, family = poisson, offset=lib.sizes)
  mu.hat.v = fitted(poi.fit)
  res.v = resid(poi.fit, type="pearson")
  
  # save as a list
  model_nb2_v_obj = list(mu.hat.v = mu.hat.v,
                         res.vec = res.v
                         )
  return(model_nb2_v_obj)
}
