
#' @title Modeling NBSTEP Genewise Dispersion with the Maximum Ajdusted Profile Likelihood (MAPL) Estimator 
#' or Maximum Likelihood Estimator (MLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NBSTEP dispersion model where the dispersion is modeled as a step (piecewise constant) function. 
#' The output of this function will be passed to the main GOF function \code{\link{nb.gof.m}}.
#' 
#' @details See the \code{\link{estimate.dispersion}} function in the \code{NBPSeq} package
#' for more information.
#' 
#' @usage
#' model.nbstep.m(counts, x, lib.sizes=colSums(counts), method=method)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param method method for estimating dispersions.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.m}} function.
#' 
#' @author Gu Mi <neo.migu@gmail.com>, Yanming Di, Daniel Schafer
#' 
#' @references 
#' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
#' Model for Assessing Differential Gene Expression from RNA-Seq", \emph{Statistical 
#' Applications in Genetics and Molecular Biology}, 10 (1).
#'  
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model.nbstep.m = function(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), method=method){
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  nc = dim(counts)[2]
  
  # preconditions
  stopifnot(is.matrix(x), nc == dim(x)[1])
  
  # data preparations
  nb.data = prepare.nb.data(counts, lib.sizes=lib.sizes)
  fit = estimate.dispersion(nb.data, x, model = "NB2", method = method, print.level=0)  # "NB2" calls disp.step() in estimate.dispersion()
  
  # extract quantities
  phi = fit$estimates        
  mu = irls.nb(y = nb.data$counts,
               s = nb.data$eff.lib.sizes, 
               x = x,
               phi = phi,
               beta0 = rep(NA, dim(x)[2]))$mu
  v = mu + phi * mu^2              # variance matrix
  res.m = as.matrix((counts - mu) / sqrt(v))  # res. matrix
  
  # make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
  res.m[ is.nan(res.m) ] = 0
  res.m[ is.infinite(res.m) ] = 0
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_nbstep_m_obj = list(mu.hat.mat = mu,
                            res.mat = res.m,
                            res.omat = res.om,
                            ord.res.vec = ord.res.v,
                            phi.hat.mat = phi
  )
  return(model_nbstep_m_obj)
}


