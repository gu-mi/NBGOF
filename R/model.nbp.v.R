
#' @title Modeling NBP Regression Model with Maximum Likelihood (ML) or Adjusted Profile
#' Likelihood (APL) on Original and Simulated Datasets
#' 
#' @description This function is designed to fit an NBP regression model. The output of
#' this function will be passed to the main GOF function \code{\link{nb.gof.v}}.
#' 
#' @details The \code{glm.nbp.1.MLE} function is used for NBP model fitting with MLE.
#' The \code{glm.nbp.1} function is used for NBP model fitting with APLE.
#' 
#' @usage
#' model.nbp.v(y, x, lib.sizes=NULL, method="ML")
#' 
#' @param y an n-by-1 vector of non-negative integers. For a typical RNA-Seq experiment, 
#' this may represent the read counts for a single gene across n samples.
#' @param x an n-by-p design matrix. If an intercept is desired in the model, you need to specify
#' the first column of \code{x} as a vector of 1.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is 1 for all samples.
#' @param method estimation method for NBP model fit. Either maximum likelihood (\code{ML})
#' or adjusted profile likelihood (\code{APL}). Default is \code{ML}.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb.gof.v}} function.
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
model.nbp.v = function(y, x, lib.sizes=NULL, method = "ML"){
  
  n = length(y)
  p = dim(x)[2]
  lib.sizes = ifelse(rep(is.null(lib.sizes),n), rep(1,n), lib.sizes)
  
  # preconditions
  stopifnot(n == dim(x)[1], n == length(lib.sizes), method %in% c("ML", "APL"))
  
  if (method == "ML"){
    # fit NBP model using MLE: glm.nbp.1.MLE()
    fit = glm.nbp.1.MLE(y=y, s=lib.sizes, x=x, print.level=0)
    mu = fit$mu
    phi = fit$phi
    res.v = fit$p.res
  }

  else if (method == "APL"){
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
