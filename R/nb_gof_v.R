
## Main function for GOF test for univariate response

################################################################################
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on
#' the regression of a univariate count response on one or more explanatory variables
#' 
#' @description This function is designed to test goodness-of-fit of an NB2 or NBP
#' negative binomial regression model with a known design matrix. 
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
#' 
#' @return An object of class "gofv" to which other methods (plot, summary, etc.)
#' can be applied.
#' 
#' @details When the response is a vector, we can use this function to test
#' the goodness-of-fit of a specified negative binomial regression model. 
#' 
#' @usage 
#' nb_gof_v(y, x, lib.sizes=NULL, sim=99, model.fit = "NB2")
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
#' seed <- 31513
#' n <- 30
#' beta.v <- c(-2,-1)
#' 
#' ## specify a design matrix
#' X <- cbind(rep(1,n),seq(1,4,length=n))
#' s <- rep(1e5,n)
#' mu <- s * exp(X %*% beta.v)
#' pi <- mu/s  # relative frequency
#' sim <- 99
#' 
#' ## simulate an n-by-1 vector for the read counts of a single gene
#' set.seed(seed)
#' alpha1 = -1  # NB1 data
#' phi0 = 0.002
#' alpha0 <- log(phi0); 
#' phi.nb1 <- phi0 * pi^alpha1
#' y.nb1 <- rnbinom(n, size=1/phi.nb1, mu=mu)
#' 
#' ## implement the GOF test for testing an NB regression model adequacy
#' # pdf("gofv-result.pdf", width=14, height=7)
#' par(mfrow=c(1,2))
#' # NB2 model fit
#' gf.nb1.nb2 <- nb_gof_v(y.nb1, X, s, sim=sim, model.fit="NB2")
#' plot(gf.nb1.nb2, conf.env=0.9, data.note = "NB1", pch=".", cex=5)
#' # NBP model fit
#' gf.nb1.nbp <- nb_gof_v(y.nb1, X, s, sim=sim, model.fit="NBP")
#' plot(gf.nb1.nbp, conf.env=0.9, data.note = "NB1", pch=".", cex=5)
#' # dev.off()
#' 
nb_gof_v <- function(y, x, lib.sizes=NULL, sim=99, model.fit = "NB2"){
  
  n = length(y)
  p = dim(x)[2]
  libs = ifelse(is.null(lib.sizes), rep(1, n), lib.sizes)
  
  # preconditions
  stopifnot(model.fit %in% c("NB2","NBP"))
  
  ## initialize simulation variables
  res.sim.mat <- matrix(0, nrow = (sim+1), ncol = n)  # residual matrix
  stat.sim <- numeric(sim)   # Pearson chi-sq test stat
  stat.sim.D = numeric(sim)  # New statistics from Dan
  stat.sim.G = numeric(sim)  # orthogonal distances
  
  #### -----------------------------------------------------------------
  if (model.fit == "NB2"){
    mnb2.0 = model_nb2_v(y=y, x=x, lib.sizes=libs)  # NB2 model fit on original data
    mu.hat.v0 = mnb2.0$mu.hat.v
    phi0 = mnb2.0$phi
    res.vec0 = mnb2.0$res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.vec.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      mnb2.h = model_nb2_v(y=y.vec.h, x=x, lib.sizes=libs)  # NB2 model fit on simulated data
      res.sim.mat[i, ] = mnb2.h$res.vec
    }
    close(pb)
    res.sim.mat[(sim+1), ] = res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "NBP"){
    mnbp.0 = model_nbp_v(y=y, x=x, lib.sizes=libs)  # NBP model fit on original data
    mu.hat.v0 = mnbp.0$mu.hat.v
    phi0 = mnbp.0$phi
    res.vec0 = mnbp.0$res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.vec.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      mnbp.h = model_nbp_v(y=y.vec.h, x=x, lib.sizes=libs)  # NBP model fit on simulated data
      res.sim.mat[i, ] = mnbp.h$res.vec
    }
    close(pb)
    res.sim.mat[(sim+1), ] = res.vec0
  }
  
  #### -----------------------------------------------------------------
  
  ## row-wise sort residual matrix res.sim.mat
  ord.res.sim.mat = t(apply(res.sim.mat, 1, sort))
  ord.typ.res.sim = apply(ord.res.sim.mat, 2, median)
  
  #### -----------------------------------------------------------------
  
  ## calcualte test statistics
  #stat0 = sum(res.vec0^2)
  stat0.D = sum( abs(ord.res.sim.mat[(sim+1), ] - ord.typ.res.sim) )
  #stat0.G = sum( sqrt((ord.res.sim.mat[(sim+1), ] - ord.typ.res.sim)^2 / 2) )
  
  for (i in 1:sim){
    #stat.sim[i] = sum(res.sim.mat[i, ]^2)
    stat.sim.D[i] = sum( abs(ord.res.sim.mat[i, ] - ord.typ.res.sim) )
    #stat.sim.G[i] = sum( sqrt((ord.res.sim.mat[i, ] - ord.typ.res.sim)^2 / 2) )
  }
  
  #### -----------------------------------------------------------------
  
  ## calculate p-values
  #   pval.P = 2 * min( (sum(stat.sim >= stat0) + 1) / (sim + 1),
  #                      1 - (sum(stat.sim >= stat0) + 1) / (sim + 1) )  # may be discarded later!
  pval.D = (sum(stat.sim.D >= stat0.D) + 1) / (sim + 1)    # ONE-SIDED!!
  #   pval.G = (sum(stat.sim.G >= stat0.G) + 1) / (sim + 1)    # ONE-SIDED!!
  #   pv.P = round(pval.P, 4)
  pv.D = round(pval.D, 6)
  #   pv.G = round(pval.G, 4)
  
  ## save as a list
  gof.obj <- list(model.fit = model.fit,
                  design.mat = x,
                  samp.size = n,
                  num.pred = p,
                  #pear.pval = pv.P,
                  new.pval = pv.D,
                  #orth.pval = pv.G,
                  res.vec0 = res.vec0,
                  ord.res.sim.mat = ord.res.sim.mat,
                  ord.typ.res.sim = ord.typ.res.sim,
                  sim = sim
  )
  class(gof.obj) <- "gofv"
  return(gof.obj)
}
