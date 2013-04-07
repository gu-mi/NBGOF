
## Main function for GOF test, in the form of a count matrix

# dispersion models expand to the following:
# NB, NBP, edgeR-common, edgeR-genewise, edgeR-tagwise, edgeR-trended
# ordered sim res matrix includes the original residual vector <--> (R+1) on denom.

# update: 2013/01/13

################################################################################
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on a 
#' RNA-Seq Dataset
#' 
#' @description This function is designed to test goodness-of-fit of different 
#' negative binomial dispersion models for a RNA-Seq dataset with a known design 
#' matrix. 
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param sim number of simulations performed.
#' @param model.fit a string of characters specifying the negative binomial
#' dispersion model used to fit the data. Currently the following dispersion models
#' are available to be checked for goodness-of-fit:
#' \itemize{
#' \item NB genewise dispersion model in the \code{\link{NBPSeq}} package (\code{NB})
#' \item NBP dispersion model in the \code{\link{NBPSeq}} package (\code{NBP})
#' \item NB common dispersion model in the \code{\link{edgeR}} package (\code{edgeR-common})
#' \item NB tagwise dispersion model in the \code{\link{edgeR}} package (\code{edgeR-tagwise})
#' \item NB trended dispersion model in the \code{\link{edgeR}} package (\code{edgeR-trended})
#' }
#' Users are recommended to specify \strong{exactly} the same characters as the
#' ones in paratheses above for each dispersion model.
#' 
#' @return An object of class "gofm" to which other methods (plot, summary, etc.)
#' can be applied.
#' 
#' @details When the response is a count matrix, we can use this function to test
#' the goodness-of-fit of a specified negative binomial dispersion model. 
#' 
#' @usage 
#' nb_gof_m(counts, x, lib.sizes=colSums(counts), sim=199, model.fit="NB")
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
#' @examples 
#' ## load package into session
#' library(NBGOF)
#' 
#' ## basic set-up of the model
#' seed = 31513
#' sim = 999
#' conf.env = 0.95
#' n = 100;
#' r = 14;
#' mu = matrix(seq(10, 1000, length=n), n, r); 
#' lib.sizes = rep(1e6, r);
#' pi = mu/lib.sizes
#' 
#' ## simulate an m-by-n count matrix mimicking a RNA-Seq dataset
#' set.seed(seed)
#' alpha1 = 0   # NB2 regression model
#' phi0 <- 0.8  # == phi
#' alpha0 = log(phi0); 
#' phi.nb2c = phi0 * pi^alpha1;  # == phi0
#' cbind(mu[,1], phi.nb2c[,1])   # make sure phi's are in reasonable range
#' y = rnbinom(n * r, mu=mu, size=1/phi.nb2c)  # response matrix
#' dim(y) = dim(mu); 
#' rownames(y) = paste("g",seq(1,n),sep="")
#' colnames(y) = paste("s",seq(1,r), sep="")
#' 
#' ## specify a design matrix
#' grp.ids = as.factor(c(rep(1,8), rep(2,4), rep(3,2)))
#' x = model.matrix(~grp.ids)
#' 
#' ## implement the GOF test for testing an NB dispersion model adequacy
#' ## CAUTION: may be time-consuming depending on the size of data and simulations
#' 
#' ## consider all dispersion estimation methods:
#' # pdf(file=file.path(path1,"gof-nb2comphi-95.pdf"), width=12, height=8)
#' par(mfrow=c(2,3))
#' #
#' fnb2.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="NB")
#' plot(fnb2.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#' #
#' fnbp.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="NBP")
#' plot(fnbp.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#' #
#' fcom.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-common")
#' plot(fcom.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#' #
#' fgen.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-genewise")
#' plot(fgen.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#' #
#' ftag.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-tagwise")
#' plot(ftag.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#' #
#' ftrd.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-trended")
#' plot(ftrd.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#' # dev.off()
#' 
nb_gof_m <- function(counts, x, lib.sizes=colSums(counts), sim=199, model.fit="NB"){
  
  stopifnot(model.fit %in% c("NB","NBP","edgeR-common","edgeR-genewise",
                             "edgeR-tagwise","edgeR-trended"))
  
  nr = dim(counts)[1]
  nc = dim(counts)[2]
  n = nr * nc
  counts.dim = paste(dim(counts)[1],"x",dim(counts)[2])
  
  ## initialize simulation variables
  ord.res.sim.mat <- matrix(0, nrow = (sim+1), ncol = n)   # ordered residual matrix
  stat.sim <- numeric(sim)   # Pearson chi-sq test stat
  #stat.sim.T <- numeric(sim) # new test stat (by Dan)
  
  #### -----------------------------------------------------------------
  if (model.fit == "NB"){
    mnb2.0 = model_nb_m(counts, x, lib.sizes=lib.sizes)
    mu.hat.mat0 = mnb2.0$mu.hat.mat
    phi.hat.mat0 = mnb2.0$phi.hat.mat
    res.omat0 = mnb2.0$res.omat
    ord.res.vec0 = mnb2.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=n, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      mnb2.h = model_nb_m(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mnb2.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "NBP"){
    mnbp.0 = model_nbp_m(counts, x, lib.sizes=lib.sizes)
    mu.hat.mat0 = mnbp.0$mu.hat.mat
    phi.hat.mat0 = mnbp.0$phi.hat.mat
    res.omat0 = mnbp.0$res.omat
    ord.res.vec0 = mnbp.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=n, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      mnbp.h = model_nbp_m(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mnbp.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-common"){
    mcom.0 = model_edgeR_common(counts, x, lib.sizes=lib.sizes)
    mu.hat.mat0 = mcom.0$mu.hat.mat
    phi.hat.mat0 = mcom.0$phi.hat.mat
    res.omat0 = mcom.0$res.omat
    ord.res.vec0 = mcom.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=n, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      mcom.h = model_edgeR_common(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mcom.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-genewise"){
    mgen.0 = model_edgeR_genewise(counts, x, lib.sizes=lib.sizes)
    mu.hat.mat0 = mgen.0$mu.hat.mat
    phi.hat.mat0 = mgen.0$phi.hat.mat
    res.omat0 = mgen.0$res.omat
    ord.res.vec0 = mgen.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=n, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      mgen.h = model_edgeR_genewise(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mgen.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-tagwise"){
    mtag.0 = model_edgeR_tagwise(counts, x, lib.sizes=lib.sizes)
    mu.hat.mat0 = mtag.0$mu.hat.mat
    phi.hat.mat0 = mtag.0$phi.hat.mat
    res.omat0 = mtag.0$res.omat
    ord.res.vec0 = mtag.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=n, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      mtag.h = model_edgeR_tagwise(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mtag.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-trended"){
    mtrd.0 = model_edgeR_trended(counts, x, lib.sizes=lib.sizes)
    mu.hat.mat0 = mtrd.0$mu.hat.mat
    phi.hat.mat0 = mtrd.0$phi.hat.mat
    res.omat0 = mtrd.0$res.omat
    ord.res.vec0 = mtrd.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=n, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      mtrd.h = model_edgeR_trended(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mtrd.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  
  typ.res.sim = apply(ord.res.sim.mat, 2, median)
  
  ## calcualte test statistics
  stat0 = sum(res.omat0^2)
  #stat0.T = sum((ord.res.vec0 - typ.res.sim)^2)
  for (i in 1:sim){
    stat.sim[i] = sum(ord.res.sim.mat[i, ]^2)
    #stat.sim.T[i] = sum((ord.res.sim.mat[i, ] - typ.res.sim)^2)
  }
  
  ## calculate p-values
  pval.P <- 2 * min( (sum(stat.sim >= stat0) + 1 ) / (sim + 1),
                        1 - ( sum(stat.sim >= stat0) + 1 ) / (sim + 1) )
  #pval.T <- 2 * min( (sum(stat.sim.T >= stat0.T) + 1 ) / (sim + 1),
  #                      1 - (sum(stat.sim.T >= stat0.T) + 1) / (sim + 1) )
  pv.P <- round(pval.P, 4)
  #pv.T <- round(pval.T, 4)
  
  ## save as a list
  gof.obj <- list(model.fit = model.fit,
                  counts.dim = counts.dim,
                  design.mat = x,
                  lib.sizes = lib.sizes,
                  pear.pval = pv.P,
                  #new.pval = pv.T,
                  #stat.sim = stat.sim,
                  #stat.sim.T = stat.sim.T,
                  mu.hat.m0 = mu.hat.mat0,
                  o.res.sim = ord.res.sim.mat,
                  typ.res.sim = typ.res.sim,
                  o.res0 = ord.res.vec0,
                  #stat0 = stat0,
                  #stat0.T = stat0.T,
                  sim = sim)
  
  # save the object as a "gofm" class
  class(gof.obj) <- "gofm"
  return(gof.obj)
}
