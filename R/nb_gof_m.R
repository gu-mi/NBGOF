## Main function for GOF test for multivariate response

# dispersion models expand to the following:
# NB, NBP, edgeR-common, edgeR-genewise, edgeR-tagwise, edgeR-trended
# ordered sim res matrix includes the original residual vector <--> (R+1) on denom.

# even though we can eliminate zero counts in original dataset, there is no guarantee that
# simulated datasets do not have zero counts. In order to proceed, we have to change those 
# genes with all zero counts into something else, say assigning the first replicate as 1

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
#' rownames(y) = paste("g", seq(1,n), sep="")
#' colnames(y) = paste("s", seq(1,r), sep="")
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
#' ## summarize the GOF test results:
#' summary(fnb2.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
#' summary(fnbp.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
#' summary(fcom.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
#' summary(fgen.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
#' summary(ftag.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
#' summary(ftrd.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")

nb_gof_m <- function(counts, x, lib.sizes=colSums(counts), sim=199, model.fit="NB"){
  
  stopifnot(model.fit %in% c("NB","NBP","edgeR-common","edgeR-genewise",
                             "edgeR-tagwise","edgeR-trended"))
  
  # we assume that genes with all zero counts have been removed to avoid zero variance
  # the first argument, counts, should be already subsetted
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  N = m * n
  counts.dim = paste(m,"x",n)
  
  ## initialize simulation variables
  ord.res.sim.mat <- matrix(0, nrow = (sim+1), ncol = N)   # ordered residual matrix
  stat.sim.G <- numeric(sim)   # orthogonal distances
  
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
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      y.mat.h[rowSums(y.mat.h) == 0, 1] = 1   # subjective!
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
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      y.mat.h[rowSums(y.mat.h) == 0, 1] = 1   # subjective!
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
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      y.mat.h[rowSums(y.mat.h) == 0, 1] = 1   # subjective!
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
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      y.mat.h[rowSums(y.mat.h) == 0, 1] = 1   # subjective!
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
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      y.mat.h[rowSums(y.mat.h) == 0, 1] = 1   # subjective!
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
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      y.mat.h[rowSums(y.mat.h) == 0, 1] = 1   # subjective!
      mtrd.h = model_edgeR_trended(y.mat.h, x, lib.sizes=lib.sizes)
      ord.res.sim.mat[i, ] = mtrd.h$ord.res.vec
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
  }
  #### -----------------------------------------------------------------
  
  # find the median of the big residual matrix (ordered): a vector
  ord.typ.res.sim = apply(ord.res.sim.mat[1:sim, ], 2, median)  # on simulated datasets ONLY!
  # subtract the typical residual vector from each row of the ordered residual matrix
  dists.mat.res =  abs(sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))
  
  # construct new distance measure matrix of dimension R-by-m
  grp.vec = ( seq_len( ncol(dists.mat.res) ) - 1 ) %/% n     # grouping vector
  dist.mat = t( rowsum(t(dists.mat.res), grp.vec) )    # orthogonal distance matrix (sim. + obs.)
  # THIS dist.mat IS UN-SORTED!! WE CAN USE THIS MATRIX FOR THE ENVELOPE METHOD CALCULATIONS!!
  
  
  # sort each row of the orthogonal distance matrix, and get typical distance vector (from sim only)
  ord.dist.mat = t(apply(dist.mat, 1, sort))
  ord.typ.dist = apply(ord.dist.mat[1:sim, ], 2, median)   # x-axis (from sim. ONLY!)
  dist.obs = ord.dist.mat[(sim+1), ]    # y-axis
  
  #### -----------------------------------------------------------------
  ## calcualte test statistics and p-values
  stat0.G = sum(ord.dist.mat[(sim+1), ])
  for (i in 1:sim){
    stat.sim.G[i] = sum(ord.dist.mat[i, ])
  }
  pval.G = (sum(stat.sim.G >= stat0.G) + 1) / (sim + 1)    # ONE-SIDED!!
  pv.G = round(pval.G, 6)
  
  #### -----------------------------------------------------------------
  ## save as a list
  gof.obj <- list(model.fit = model.fit,
                  counts.dim = counts.dim,
                  design.mat = x,
                  lib.sizes = lib.sizes,
                  orth.pval = pv.G,
                  mu.hat.m0 = mu.hat.mat0,
                  ord.dist.mat = ord.dist.mat,
                  ord.typ.dist  = ord.typ.dist,
                  dist.obs = dist.obs,
                  dist.mat = dist.mat,
                  sim = sim)
  
  # save the object as a "gofm" class
  class(gof.obj) <- "gofm"
  return(gof.obj)
}
