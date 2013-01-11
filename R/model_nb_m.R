
## For NB2 model fitting on the original & simulated datasets
## use nb.regression.1() function for NB2 fitting 
## (avoid convergence issues in MASS::glm.nb)

# update: 2013/01/08

#' @title Modeling NB genewise model with MLE on original and simulated datasets
#' 
#' @description This function is designed to fit an NB regression model with
#' genewise dispersions using the maximum likelihood estimator. The output of
#' this function will be passed to the main GOF function
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples
#' @param x an n-by-p design matrix
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix
#' 
#' @return An object called "model_nb_m_obj" to be passed to the main
#' \code{\link{nb_gof_m}} function
#' 
#' @author Gu Mi <\url{http://people.oregonstate.edu/~mig}>
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_nb_m <- function(counts, x, lib.sizes = colSums(counts)){
  
  nc = dim(counts)[2]  
  
  # preconditions
  stopifnot(is.matrix(x), nc == dim(x)[1])
  
  # NB2 fit and extract quantities
  nb.fit = nb.regression.1(y=counts, s=lib.sizes, x=x, beta=NA)
  mu.hat.m = nb.fit$mu
  v.hat.m = nb.fit$v
  phi.hat.m = nb.fit$phi
  res.m = (counts - mu.hat.m) / sqrt(v.hat.m) 
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_nb_m_obj = list(mu.hat.mat = mu.hat.m,
                        res.mat = res.m,
                        res.omat = res.om,
                        ord.res.vec = ord.res.v,
                        phi.hat.mat = phi.hat.m,
                        fit.nb = nb.fit
  )
  return(model_nb_m_obj)
}
