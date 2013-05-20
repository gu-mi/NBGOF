
## Main function for GOF test for univariate response

################################################################################
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on
#' the regression of a univariate count response on one or more explanatory variables
#' 
#' @description This function is designed to test goodness-of-fit of an NB2 or NBP
#' negative binomial regression model with a known design matrix. Estimation methods
#' for NBP model fitting include MLE and APLE
#' 
#' @param y an n-by-1 vector of non-negative integers. For a typical RNA-Seq experiment, 
#' this may represent the read counts for a single gene 
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is 1 for all samples
#' @param sim number of simulations performed.
#' @param model.fit a string of characters specifying the negative binomial model used 
#' to fit the data. Currently the following dispersion models are available to be checked
#' for goodness-of-fit:
#' \itemize{
#' \item NB2: conventional NB2 regression model (\code{NB2})
#' \item NBP: NBP regression model (see \code{\link{NBPSeq}}) package (\code{NBP})
#' }
#' Users are recommended to specify \strong{exactly} the same characters as the
#' ones in paratheses above for each NB regression model.
#' @param est.method either "MLE" or "APLE" for the estimation methods used
#' 
#' @return An object of class "gofv" to which other methods (plot, summary, etc.)
#' can be applied.
#' 
#' @details When the response is a vector, we can use this function to test
#' the goodness-of-fit of a specified negative binomial regression model. 
#' 
#' @usage 
#' nb_gof_v(y, x, lib.sizes=NULL, sim=999, model.fit = "NB2", est.method="MLE")
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
#' @examples 
#' 
#' ## load package into session
#' library(NBGOF)
#' 
#' ## basic set-up of the model
#' seed = 31513
#' n = 100
#' beta.v = c(1, -3)
#' 
#' ## specify a design matrix
#' X = cbind(rep(1,n),seq(1.5, 3.5, length=n))
#' s = rep(1e6, n)
#' mu = s * exp(X %*% beta.v)
#' pi = mu/s  # relative frequency
#' sim = 499
#' 
#' ## simulate an n-by-1 vector for the read counts of a single gene
#' set.seed(seed)
#' alpha1 = -1  # NB1 data
#' phi0 = 0.0001
#' alpha0 = log(phi0); 
#' phi.nb1 = phi0 * pi^alpha1
#' y.nb1 = rnbinom(n, size=1/phi.nb1, mu=mu)
#' 
#' ## implement the GOF test for testing an NB regression model adequacy of NB2 and NBP models
#' # pdf("gofv-result.pdf", width=14, height=7)
#' par(mfrow=c(1,2))
#' # NB2 model fit
#' gf.nb1.nb2 = nb_gof_v(y.nb1, X, s, sim=sim, model.fit="NB2")
#' plot(gf.nb1.nb2, conf.env=0.95, data.note = "NB1", pch=".", cex=5)
#' # NBP model fit
#' gf.nb1.nbp = nb_gof_v(y.nb1, X, s, sim=sim, model.fit="NBP")
#' plot(gf.nb1.nbp, conf.env=0.95, data.note = "NB1", pch=".", cex=5)
#' # dev.off()
#' 
nb_gof_v = function(y, x, lib.sizes=NULL, sim=999, model.fit = "NB2", est.method="MLE"){
  
  n = length(y)
  p = dim(x)[2]
  
  # preconditions
  stopifnot(model.fit %in% c("Poisson", "NB2", "NBP"))
  
  ## initialize simulation variables
  res.sim.mat = matrix(0, nrow = (sim+1), ncol = n)  # residual matrix
  stat.sim.P = numeric(sim)  # Pearson chi-sq test stat (NOT USED! for simulation only!)
  stat.sim.D = numeric(sim)  # statistic based on the overall vertical distance
  
  #### -----------------------------------------------------------------
  if (model.fit == "Poisson"){
    libs = ifelse(is.null(lib.sizes), rep(0, n), lib.sizes)  # pay attention to the offset here!
    mpoi.0 = model_poi_v(y=y, x=x, lib.sizes=libs)  
    # Poisson model fit on original data
    mu.hat.v0 = mpoi.0$mu.hat.v
    res.vec0 = mpoi.0$res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.vec.h = rpois(n, lambda=mu.hat.v0)
      mpoi.h = model_poi_v(y=y.vec.h, x=x, lib.sizes=libs)  
      # Poisson model fit on simulated data
      res.sim.mat[i, ] = mpoi.h$res.vec
    }
    close(pb)
    res.sim.mat[(sim+1), ] = res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "NB2"){
    libs = ifelse(is.null(lib.sizes), rep(1, n), lib.sizes)
    mnb2.0 = model_nb2_v(y=y, x=x, lib.sizes=libs)  
    # NB2 model fit on original data
    mu.hat.v0 = mnb2.0$mu.hat.v
    phi0 = mnb2.0$phi
    res.vec0 = mnb2.0$res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.vec.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      mnb2.h = model_nb2_v(y=y.vec.h, x=x, lib.sizes=libs)  
      # NB2 model fit on simulated data
      res.sim.mat[i, ] = mnb2.h$res.vec
    }
    close(pb)
    res.sim.mat[(sim+1), ] = res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "NBP"){
    libs = ifelse(is.null(lib.sizes), rep(1, n), lib.sizes)
    mnbp.0 = model_nbp_v(y=y, x=x, lib.sizes=libs, est.method=est.method)  
    # NBP model fit on original data
    mu.hat.v0 = mnbp.0$mu.hat.v
    phi0 = mnbp.0$phi
    res.vec0 = mnbp.0$res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.vec.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      mnbp.h = model_nbp_v(y=y.vec.h, x=x, lib.sizes=libs, est.method=est.method)  
      # NBP model fit on simulated data
      res.sim.mat[i, ] = mnbp.h$res.vec
    }
    close(pb)
    res.sim.mat[(sim+1), ] = res.vec0
  }
  
  #### -----------------------------------------------------------------
  ## row-wise sort residual matrix res.sim.mat
  ord.res.sim.mat = t(apply(res.sim.mat, 1, sort))
  ord.typ.res.sim = apply(ord.res.sim.mat[-(sim+1), ], 2, median)   # exclude the obs. row
  
  #### -----------------------------------------------------------------
  ## calcualte test statistics (observed and simulated) and p-values
  stat0.P = sum(res.vec0^2)
  stat0.D = sum( (ord.res.sim.mat[(sim+1), ] - ord.typ.res.sim)^2 )
  
  for (i in 1:sim){
    stat.sim.P[i] = sum( res.sim.mat[i, ]^2 )
    stat.sim.D[i] = sum( (ord.res.sim.mat[i, ] - ord.typ.res.sim)^2 )
  }
  pval.P = (sum(stat.sim.P >= stat0.P) + 1) / (sim + 1)
#   pval.P = 2 * min( (sum(stat.sim.P >= stat0.P) + 1 ) / (sim + 1),
#                     1 - ( sum(stat.sim.P >= stat0.P) + 1 ) / (sim + 1) )  # TWO-SIDED
  pval.D = (sum(stat.sim.D >= stat0.D) + 1) / (sim + 1)
  pv.P = round(pval.P, 6)
  pv.D = round(pval.D, 6)
  
  #### -----------------------------------------------------------------
  ## save as a list
  gof.obj = list(model.fit = model.fit,
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
  class(gof.obj) = "gofv"
  return(gof.obj)
}
