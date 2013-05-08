
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on
#' the regression of a univariate count response on one or more explanatory variables
#' 
#' @description This function is designed to test goodness-of-fit of an NB2 or NBP
#' negative binomial regression model with a known design matrix. Estimation methods
#' for NBP model fitting include MLE and APLE
#' 
#' @param y0 an n-by-1 vector of non-negative integers. For a typical RNA-Seq experiment, 
#' this may represent the read counts for a single gene 
#' @param x an n-by-p design matrix.
#' @param sim number of simulations performed.
#' @param model a string of characters specifying the negative binomial model used 
#' to fit the data. Currently the following dispersion models are available to be checked
#' for goodness-of-fit
#' 
#' @return An object
#' 
#' @details When the response is a vector, we can use this function to test
#' the goodness-of-fit of a specified negative binomial regression model. 
#' 
#' @usage 
#' ks_gof(y0, x, model = "NB2", sim=1)
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 

ks_gof = function(y0, x, model = "NB2", sim=1){
  
  if (model == "NB2"){
    
    n = length(y0)
    res.h.mat = matrix(0, nrow = sim, ncol = n)
    
    # get Pearson residuals on observed response y0:
    fit.nb2 = model_nb2_v(y=y0, x=x)
    res0 = fit.nb2$res.vec
    # ordered Pearson residuals:
    ord.res0 = sort(res0)
    
    # extract fitted quantities:
    mu.hat.v0 = fit.nb2$mu.hat.v
    phi0 = fit.nb2$phi
    
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      # re-simulate new responses y:
      y.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      fit.nb2.h = model_nb2_v(y=y.h, x=x)
      res.h.mat[i, ] = fit.nb2.h$res.vec
    }
    close(pb)
    
    res.sort = t(apply(res.h.mat, 1, sort))
    ord.res.h = apply(res.sort, 2, median)
    #
    ks_test_res = ks.test(x=ord.res0, y=ord.res.h, alternative="two.sided")  
  }
  return(list(pval = ks_test_res$p.value,
              ord.res0 = ord.res0,
              ord.res.h = ord.res.h)
  )
  
}