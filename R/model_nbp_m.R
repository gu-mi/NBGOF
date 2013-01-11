
## For NBP model fitting on the original & simulated dataset

# update: 2013/01/08

#' @title Modeling NBP genewise model with the maximum likelihood estimator (MLE) on 
#' original and simulated datasets
#' 
#' @description This function is designed to fit an NBP regression model with
#' genewise dispersions using the adjusted profile likelihood estimator. See details
#' below. The output of this function will be passed to the main GOF function
#' 
#' @details Under the NB model, the mean-variance relationship of a single read count 
#' satisfies $\sigma_{ij}^2 = \mu_{ij} + \phi_{ij} \mu_{ij}^2$. For applying the NBP 
#' model to RNA-Seq data, we consider the "log-linear-rel-mean" method assuming a 
#' parametric dispersion model
#' \[
#' \ohi_{ij} = \alpha_0 + \alpha_1 \log(\pi_{ij}),
#' \]
#' where $\pi_{ij} = \mu_{ij}/(N_j R_j) is the relative mean frequency after 
#' normalization. The parameters $(\alpha_0, \alpha_1)$ in this dispersion model 
#' are estimated by maximizing the adjusted profile likelihood. See the 
#' \code{\link{estimate.dispersion}} function in the \code{\link{NBPSeq}} package
#' for more information.
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples
#' @param x an n-by-p design matrix
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix
#' 
#' @return An object called "model_nbp_m_obj" to be passed to the main
#' \code{\link{nb_gof_m}} function
#' 
#' @author Gu Mi <\url{http://people.oregonstate.edu/~mig}>
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
model_nbp_m <- function(counts, x, lib.sizes = colSums(counts)){
  
  nc = dim(counts)[2]
  
  # preconditions
  stopifnot(is.matrix(x), nc == dim(x)[1])
  
  # data preparations
  nb.data <- prepare.nb.data(counts, lib.sizes=lib.sizes)
  fit <- estimate.dispersion(nb.data, x, print.level=0)
  phi <- fit$models[[1]]$phi        # NBP "phi" --> mu+phi*mu^2
  mu <- fit$models[[1]]$mu
  v <- mu + phi * mu^2              # variance matrix
  res.m <- (counts - mu) / sqrt(v)  # res. matrix
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_nbp_m_obj = list(mu.hat.mat = mu,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi,
                         fit.nbp = fit$models[[1]]
  )
  return(model_nbp_m_obj)
}

