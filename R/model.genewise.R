
#' @title Modeling NB Genewise Dispersion Using NBPSeq
#' 
#' @description This function fits an NB regression model with
#' genewise dispersions using the adjusted profile likelihood estimator. See details below. The output of this function 
#' will be passed to the main GOF function \code{\link{nb.gof.m}}.
#' 
#' @details details here (HOA)
#' 
#' @usage
#' model.genewise(counts, x)
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.m}} function.
#'  
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
model.genewise = function(counts, x){
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  
  # keep only complete cases:
  counts = counts[complete.cases(counts), ]
  
  nf = estimate.norm.factors(counts, method="AH2010")
  nb.data = prepare.nb.data(counts, norm.factors=nf)
  grp1 = as.character(unique(grp.ids)[1])
  grp2 = as.character(unique(grp.ids)[2])
  gen.fit = genewise(nb.data=nb.data, grp.ids=grp.ids, grp1=grp1, grp2=grp2, R = 100)
  
  # extract quantities:
  mu.hat.m = matrix(gen.fit$mu.hat, ncol = length(grp.ids))   # mu may be close to 0
  phi.hat.m = gen.fit$phi.hat     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2
  res.m = as.matrix((counts - mu.hat.m) / sqrt(v))
  
  # make sure 0/0 (NaN) and 1/0 (Inf) won't appear in residual matrix (before sorting)
  res.m[ is.nan(res.m) ] = 0
  res.m[ is.infinite(res.m) ] = 0
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort.vec, grp.ids)) 
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_gen_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
  )
  return(model_gen_m_obj)
}




