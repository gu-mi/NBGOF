
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on
#' the Regression of a Univariate Count Response on One or More Explanatory Variables
#' 
#' @description This function tests goodness-of-fit of an NB2 or NBP
#' negative binomial regression model with a known design matrix. Estimation method 
#' for the NB2 model fitting is maximum likelihood (ML). Estimation methods for NBP model fitting
#' include both ML and adjusted profile likelihood (APL). This function
#' can also test goodness-of-fit of a Poisson regression model.
#' 
#' @param y an n-by-1 vector of non-negative integers. For a typical RNA-Seq experiment, 
#' this may represent the read counts for a single gene across n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is 1 for all samples.
#' @param sim number of simulations performed.
#' @param model a string of characters specifying the model (negative binomial or Poisson) 
#' used to fit the data. Currently the following models are available to be checked
#' for goodness-of-fit:
#' \itemize{
#' \item NB2: conventional NB2 regression model (\code{NB2})
#' \item NBP: NBP regression model (\code{NBP})
#' \item Poisson: Poisson regression model (\code{Poisson})
#' }
#' Users should specify \strong{exactly} the same characters as shown in paratheses above 
#' for testing one of the regression models.
#' @param method specify either "\code{ML}" or "\code{APL}" for the maximum likelihood
#' and adjusted profile likelihood methods used, respectively, for the NBP model estimations. ML is used
#' for the NB2 model estimations.
#' 
#' @return An object of class "gofv" to which other methods can be applied.
#' 
#' @details When the response is a vector of counts, we can use this function to test
#' the goodness-of-fit of a specified negative binomial or Poisson regression model. It returns
#' an object with the test results, which can be further summarized and visualized using 
#' appropriate methods, e.g. \code{\link{EPPlot}}.
#' 
#' This function calls \code{\link{model.nb2.v}} to fit the NB2 model, calls 
#' \code{\link{model.nbp.v}} to fit the NBP model, and calls \code{\link{model.poi.v}} to fit
#' the Poisson model.
#' 
#' @usage 
#' nb.gof.v(y, x, lib.sizes=NULL, sim=999, model = "NB2", method="ML")
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references 
#' Mi G., Di Y. and Schafer D.W. (2013). Goodness-of-Fit Tests and Model Diagnostics
#' for Negative Binomial Regression of RNA sequencing Data. \emph{Biometrics} (revision invited).
#' 
#' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
#' Model for Assessing Differential Gene Expression from RNA-Seq", \emph{Statistical 
#' Applications in Genetics and Molecular Biology}, 10 (1).
#' 
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
#' @examples 
#' 
#' ## load package into session:
#' library(NBGOF)
#' 
#' ## basic set-up of the model:
#' seed = 539768
#' n = 100
#' beta.v = c(1, -3)
#' 
#' ## specify a design matrix:
#' X = cbind(rep(1,n), seq(1.5, 3.5, length=n))
#' s = rep(1e6, n)
#' mu = s * exp(X %*% beta.v)
#' pi = mu/s  # relative frequency
#' sim = 999
#' 
#' ## simulate an n-by-1 vector for the read counts of a single gene:
#' # the data is simulated from NB1
#' set.seed(seed)
#' alpha1 = -1  # NB1 data
#' phi0 = 0.0001
#' alpha0 = log(phi0); 
#' phi.nb1 = phi0 * pi^alpha1
#' y.nb1 = rnbinom(n, size=1/phi.nb1, mu=mu)
#' 
#' ## implement the GOF test for testing NB and Poisson model adequacy:
#' # pdf("gofv-result.pdf", width=14, height=7)
#' par(mfrow=c(1,3))
#' 
#' # NB2 model fit using MLE:
#' gof.nb1.nb2 = nb.gof.v(y.nb1, X, s, sim=sim, model="NB2")
#' EPPlot(gof.nb1.nb2, conf.env=0.95, data.note="NB1", pch=".", cex=5)
#' 
#' # NBP model fit using MLE:
#' gof.nb1.nbp = nb.gof.v(y.nb1, X, s, sim=sim, model="NBP", method="ML")
#' EPPlot(gof.nb1.nbp, conf.env=0.95, data.note="NB1", pch=".", cex=5)
#' 
#' # Poisson model fit:
#' gof.nb1.poi = nb.gof.v(y.nb1, X, s, sim=sim, model="Poisson")
#' EPPlot(gof.nb1.poi, conf.env=0.95, data.note="NB1", pch=".", cex=5)
#' # dev.off()
#' 
nb.gof.v = function(y, x, lib.sizes=NULL, sim=999, model = "NB2", method="ML", ncores = NULL){
  
  n = length(y)
  p = dim(x)[2]
  
  # preconditions
  stopifnot(model %in% c("Poisson", "NB2", "NBP"), method %in% c("ML","APL"))
  
  # parallel computing: specify number of cores to use
  if (is.null(ncores)){
    ncores = detectCores() - 1
  }
  # register for foreach
  registerDoMC(ncores)
  
  ## initialize simulation variables
  # res.sim.mat = matrix(0, nrow = (sim+1), ncol = n)  # residual matrix
  stat.sim.P = numeric(sim)  # Pearson chi-sq test statistic
  stat.sim.D = numeric(sim)  # statistic based on the overall vertical distance
  
  #### -----------------------------------------------------------------
  if (model == "Poisson"){
    libs = ifelse(is.null(lib.sizes), rep(0, n), lib.sizes)  # pay attention to the offset here!
    mpoi.0 = model.poi.v(y=y, x=x, lib.sizes=libs)  
    # Poisson model fit on original data
    mu.hat.v0 = mpoi.0$mu.hat.v
    res.vec0 = mpoi.0$res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ## ---------------------------------
    ## Parallel computing begins here ##
    ## ---------------------------------    
    res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      y.vec.h = rpois(n, lambda=mu.hat.v0)
      mpoi.h = model.poi.v(y=y.vec.h, x=x, lib.sizes=libs)  
      # Poisson model fit on simulated data
      mpoi.h$res.vec
    }
    #close(pb)
    dimnames(res.sim.mat.tmp) = NULL
    res.sim.mat = rbind(res.sim.mat.tmp, res.vec0)
  }
  #### -----------------------------------------------------------------
  if (model == "NB2"){
    libs = ifelse(is.null(lib.sizes), rep(1, n), lib.sizes)
    mnb2.0 = model.nb2.v(y=y, x=x, lib.sizes=libs)  # MLE for NB2 models
    # NB2 model fit on original data
    mu.hat.v0 = mnb2.0$mu.hat.v
    phi0 = mnb2.0$phi
    res.vec0 = mnb2.0$res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ## ---------------------------------
    ## Parallel computing begins here ##
    ## ---------------------------------    
    res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      y.vec.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      mnb2.h = model.nb2.v(y=y.vec.h, x=x, lib.sizes=libs)  # MLE for NB2 models
      # NB2 model fit on simulated data
      mnb2.h$res.vec
    }
    #close(pb)
    dimnames(res.sim.mat.tmp) = NULL
    res.sim.mat = rbind(res.sim.mat.tmp, res.vec0)
  }
  #### -----------------------------------------------------------------
  if (model == "NBP"){
    libs = ifelse(is.null(lib.sizes), rep(1, n), lib.sizes)
    mnbp.0 = model.nbp.v(y=y, x=x, lib.sizes=libs, method=method)  
    # NBP model fit on original data
    mu.hat.v0 = mnbp.0$mu.hat.v
    phi0 = mnbp.0$phi
    res.vec0 = mnbp.0$res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ## ---------------------------------
    ## Parallel computing begins here ##
    ## ---------------------------------    
    res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      y.vec.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      mnbp.h = model.nbp.v(y=y.vec.h, x=x, lib.sizes=libs, method=method)  
      # NBP model fit on simulated data
      mnbp.h$res.vec
    }
    #close(pb)
    dimnames(res.sim.mat.tmp) = NULL
    res.sim.mat = rbind(res.sim.mat.tmp, res.vec0)
  }
  
  #### -----------------------------------------------------------------
  ## row-wise sort residual matrix res.sim.mat ("ordered")
  ord.res.sim.mat = t(apply(res.sim.mat, 1, sort))
  ## get the typical residual vector from column medians
  ## note: we do NOT include the observed residual vector in the calculation!
  ord.typ.res.sim = apply(ord.res.sim.mat[-(sim+1), ], 2, median)
  
  #### -----------------------------------------------------------------
  ## calcualte test statistics (observed and simulated) and p-values
  stat0.P = sum(res.vec0^2)
  stat0.D = sum( (ord.res.sim.mat[(sim+1), ] - ord.typ.res.sim)^2 )
  #
  for (i in 1:sim){
    stat.sim.P[i] = sum( res.sim.mat[i, ]^2 )
    stat.sim.D[i] = sum( (ord.res.sim.mat[i, ] - ord.typ.res.sim)^2 )
  }
  # 2-sided p-value for Pearson statistics
  # 1-sided p-value for vertical distance measure
  pval.P = 2 * min( (sum(stat.sim.P >= stat0.P) + 1 ) / (sim + 1),
                    1 - ( sum(stat.sim.P >= stat0.P) + 1 ) / (sim + 1) )
  pval.D = (sum(stat.sim.D >= stat0.D) + 1) / (sim + 1)
  pv.P = round(pval.P, 6)
  pv.D = round(pval.D, 6)
  
  #### -----------------------------------------------------------------
  ## save as a list
  gof.obj = list(model = model,
                 design.mat = x,
                 samp.size = n,
                 num.pred = p,
                 pear.pval = pv.P,
                 new.pval = pv.D,
                 res.vec0 = res.vec0,
                 ord.res.sim.mat = ord.res.sim.mat,
                 ord.typ.res.sim = ord.typ.res.sim,
                 sim = sim
  )
  # save the object as a "gofv" class
  class(gof.obj) = "gofv"
  return(gof.obj)
}
