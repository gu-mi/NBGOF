
## Main function for GOF test, in the form of a count matrix

# dispersion models expand to the following:
# NB2, NBP, edgeR-common, edgeR-genewise, edgeR-tagwise, edgeR-trended
# ordered sim res matrix includes the original residual vector <--> (R+1) on denom.

# update: 2013/01/08

# library(MASS)
# library(NBPSeq)
# library(edgeR)
# library(calibrate)
# 
# source("nb.regression.1.R") # new code for NB fitting (MLE)
# source("model_nb_m.R")     # for modeling NB on original & simulated datasets
# source("model_nbp_m.R")     # for modeling NBP on original & simulated datasets
# source("model_edgeR_common.R")   # for modeling NB2 common disp. model by edgeR
# source("model_edgeR_genewise.R") # for modeling NB2 genewise disp. model by edgeR
# source("model_edgeR_tagwise.R")  # for modeling NB2 tagwise disp. model by edgeR
# source("model_edgeR_trended.R")  # for modeling NB2 tagwise disp. model by edgeR
# source("plot.gofm.R")        # for plotting an object of class "gofm"

#' @title Implement simulation-based goodness-of-fit tests on a RNA-Seq dataset
#' 
#' @description This function is designed to test goodness-of-fit of different 
#' negative binomial dispersion models for a RNA-Seq dataset with a known design 
#' matrix. 
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples
#' @param x an n-by-p design matrix
#' @param lib.sizes library sizes of a RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix
#' @param sim number of simulations performed
#' @param model.fit a string of characters specifying the negative binomial
#' dispersion model used to fit the data.
#' 
#' @return An object called "gofm" to which other methods (plot, summary, etc.)
#' can be applied.
#' 
#' @author Gu Mi <\url{http://people.oregonstate.edu/~mig}>
#' 
#' @export
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
#' @examples library(NBGOF)
#' ## simulate an m-by-n count matrix mimicking a RNA-Seq dataset
#' 
#' ## specify a design matrix
#' 
#' ## implement the GOF test for testing an NB dispersion model adequacy
#' 
nb_gof_m <- function(counts, x, lib.sizes = colSums(counts), sim=199, model.fit = "NB"){
  
  stopifnot(model.fit %in% c("NB","NBP","edgeR-common","edgeR-genewise",
                             "edgeR-tagwise","edgeR-trended"))
  
  nr = dim(counts)[1]
  nc = dim(counts)[2]
  n = nr * nc
  counts.dim = paste(dim(counts)[1],"x",dim(counts)[2])
  
  ## initialize simulation variables
  ord.res.sim.mat <- matrix(0, nr = (sim+1), nc = n)   # ordered residual matrix
  stat.sim <- numeric(sim)   # Pearson chi-sq test stat
  stat.sim.T <- numeric(sim) # new test stat (by Dan)
  
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
      mnb2.h = model_nb2_m(y.mat.h, x, lib.sizes=lib.sizes)
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
  stat0.T = sum((ord.res.vec0 - typ.res.sim)^2)
  for (i in 1:sim){
    stat.sim[i] = sum(ord.res.sim.mat[i, ]^2)
    stat.sim.T[i] = sum((ord.res.sim.mat[i, ] - typ.res.sim)^2)
  }
  
  ## calculate p-values
  pval.P <- 2 * min( (sum(stat.sim >= stat0) + 1 ) / (sim + 1),
                        1 - ( sum(stat.sim >= stat0) + 1 ) / (sim + 1) )
  pval.T <- 2 * min( (sum(stat.sim.T >= stat0.T) + 1 ) / (sim + 1),
                        1 - (sum(stat.sim.T >= stat0.T) + 1) / (sim + 1) )
  pv.P <- round(pval.P, 4)
  pv.T <- round(pval.T, 4)
  
  ## save as a list
  gof.obj <- list(model.fit = model.fit,
                  counts.dim = counts.dim,
                  design.mat = x,
                  lib.sizes = lib.sizes,
                  pear.pval = pv.P,
                  new.pval = pv.T,
                  stat.sim = stat.sim,
                  stat.sim.T = stat.sim.T,
                  mu.hat.m0 = mu.hat.mat0,
                  o.res.sim = ord.res.sim.mat,
                  typ.res.sim = typ.res.sim,
                  o.res0 = ord.res.vec0,
                  stat0 = stat0,
                  stat0.T = stat0.T,
                  sim = sim)
  
  # save the object as a "gofm" class
  class(gof.obj) <- "gofm"
  return(gof.obj)
}
