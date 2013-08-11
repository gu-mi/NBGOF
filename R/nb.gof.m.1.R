
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on an 
#' RNA-Seq Dataset
#' 
#' @description This function performs simulation-based goodness-of-fit tests of different 
#' negative binomial dispersion models for an RNA-Seq dataset with a known design 
#' matrix. Currently supported models are implemented in the \code{NBPSeq} and \code{edgeR}
#' packages.
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' @param sim number of simulations performed.
#' @param model.fit a string of characters specifying the negative binomial
#' dispersion model used to fit the data. Currently the following dispersion models
#' are available to be checked for goodness-of-fit:
#' \itemize{
#' \item NBP dispersion model in the \code{NBPSeq} package (\code{NBP})
#' \item NB common dispersion model in the \code{\link{edgeR}} package (\code{edgeR-common})
#' \item NB genewise dispersion model in the \code{\link{edgeR}} package (\code{edgeR-genewise})
#' \item NB trended (non-parametric) dispersion model in the \code{\link{edgeR}} package
#' (\code{edgeR-trended})
#' \item NB tagwise-common dispersion model in the \code{\link{edgeR}} package 
#' (\code{edgeR-tagcom})
#' \item NB tagwise-trended dispersion model in the \code{\link{edgeR}} package 
#' (\code{edgeR-tagtrd})
#' }
#' Users should specify \strong{exactly} the same characters as the
#' ones in paratheses above for each dispersion model.
#' @param min.n for \code{edgeR-trended} model only: specify the minimim number of genes in a bin
#' (default is 10).
#' @param prior.df control of the shrinkage applied for \code{edgeR} genewise, tagwise 
#' (and its variants) dispersion models. Setting \code{prior.df=0} gives the genewise model
#' without any shrinkage. Setting \code{prior.df=10} (default in \code{edgeR}) gives the tagwise
#' model with shrinkage, either towards the common dispersion (\code{edgeR-tagcom}) or the
#' global dispersion trend (\code{edgeR-tagtrd}).
#' @param design specifications of testing (1) a single-group model (\code{single}); (2)
#' a multiple-group model (\code{multiple}); (3) complex design with interactions
#' (\code{complex}). \strong{The codes are only tested for single- and multiple-group cases, and
#' the modeling of dispersions in the multiple-group case is still under consideration!
#' For the complex design with interactions, though we can pass the design matrix directly
#' in edgeR, the results may not make any sense!}
#' 
#' @return An object of class "gofm" to which other methods can be applied.
#' 
#' @details When the response is a count matrix, we can use this function to test
#' the goodness-of-fit of a specified negative binomial dispersion model. Specifically,
#' this function calls \code{\link{model_nbp_m}} to fit the NBP dispersion model, calls
#' \code{\link{model_edgeR_common}} to fit the common dispersion model, calls
#' \code{\link{model_edgeR_genewise}} to fit the pure genewise dispersion model (without
#' shrinkage), calls \code{\link{model_edgeR_trended}} to fit the non-parametric trended
#' dispersion model in \code{edgeR}, calls \code{\link{model_edgeR_tagcom}} to fit the tagwise
#' model with the dispersions shrinking towards a common dispersion, and calls
#' \code{\link{model_edgeR_tagtrd}} to fit the tagwise model with the dispersions shrinking
#' towards the global dispersion trend.
#' 
#' @usage 
#' nb_gof_m(counts, x, lib.sizes=colSums(counts), sim=999, model.fit="NBP", min.n=100, 
#' prior.df = 10, design = "simple")
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references 
#' Mi G, Di Y and Schafer DW (2013). Goodness-of-Fit Tests and Model Diagnostics
#' for Negative Binomial Regression of RNA sequencing Data. \emph{Biometrics} (submitted).
#' 
#' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
#' Model for Assessing Differential Gene Expression from RNA-Seq", \emph{Statistical 
#' Applications in Genetics and Molecular Biology}, 10 (1).
#' 
#' McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor 
#' RNA-Seq experiments with respect to biological variation. \emph{Nucleic Acids Research} 
#' 40, 4288-4297.
#' 
#' Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference.
#' \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
#' 
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
#' @examples 
#' ## load package into session:
#' library(NBGOF)
#' 
#' ## basic set-up of the model:
#' seed = 315825
#' sim = 999
#' conf.env = 0.95
#' 
#' ## basic parameter specifications:
#' m = 1000; # 1000 genes
#' n = 6;    # 6 samples
#' s = 1e6;  # lib.sizes
#' offset = log(s);  # make sure offset for NBP functions should take the log
#' offset    # 13.82
#' 
#' ## simulate coefficients:
#' beta = matrix(0, m, 1);  # beta0 only
#' set.seed(seed)
#' beta[,1] = rnorm(m, 5, 2) - offset;   # beta0 ~ normal (mean and a based on real RNA-Seq data)
#' beta  # m-by-1 matrix
#' 
#' ## design matrix (single group only):
#' x = matrix(rep(1,n))
#' x
#' 
#' ## specify mean levels:
#' mu = round(t(s * exp(x %*% t(beta))));
#' sum(rowSums(mu) == 0)   # 4
#' mu[rowSums(mu) == 0, ] = 1
#' sum(rowSums(mu) == 0)   # 0
#' pi = mu/s
#' 
#' ## simulate an m-by-n count matrix mimicking a RNA-Seq dataset:
#' ## simulate NB1.8 data with random uniform noise added to the dispersion
#' alpha1 = -0.2  # NB1.8
#' phi0 = 0.02
#' alpha0 = log(phi0)
#' phi.nbp = phi0 * pi^alpha1
#' range(phi.nbp)
#' cbind(mu[,1], phi.nbp[,1]);
#' 
#' ## add noise:
#' a = 0.1
#' set.seed(seed)  
#' phi.noi.vec = phi.nbp[ ,1] * exp(matrix(runif(m, -a, a), nr=m, nc=1))
#' phi.noi = matrix(phi.noi.vec, nr=m, nc=n)
#' range(phi.noi)
#' cbind(mu[,1], phi.noi[,1])  # make sure phi's are in reasonable range
#' 
#' ## generate NBP response with added noise:
#' set.seed(seed)
#' y = rnbinom(m * n, mu=mu, size=1/phi.noi)
#' dim(y) = dim(mu)
#' rownames(y) = paste("g", seq(1,m), sep="")
#' colnames(y) = paste("s", seq(1,n), sep="")
#' plot(mu, phi.noi, log="xy")
#' 
#' ## GOF tests for different dispersion models:
#' ## CAUTION: may be time-consuming depending on the size of data and simulations
#' fnbp.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="NBP")
#' fcom.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-common")
#' fgen.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-genewise")
#' ftgc.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-tagcom")
#' ftgt.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-tagtrd")
#' ftrd.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-trended")
#' 
#' ## summarize the GOF test results:
#' summary(fnbp.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(fcom.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(fgen.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(ftgc.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(ftgt.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(ftrd.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' 
#' ## evaluate GOF p-values using histograms and quantile-quantile uniform plots:
#' # pdf(file="plots-NBP-01-1grp.pdf", width=8, height=24)
#' par(mfrow=c(6,2))
#' plot(fnbp.noip, model="NBP", type="hist")
#' plot(fnbp.noip, model="NBP", type="qq")
#' plot(fcom.noip, model="Common", type="hist")
#' plot(fcom.noip, model="Common", type="qq")
#' plot(fgen.noip, model="Genewise", type="hist")
#' plot(fgen.noip, model="Genewise", type="qq")
#' plot(ftgc.noip, model="Tagwise-Common", type="hist")
#' plot(ftgc.noip, model="Tagwise-Common", type="qq")
#' plot(ftgt.noip, model="Tagwise-Trended", type="hist")
#' plot(ftgt.noip, model="Tagwise-Trended", type="qq")
#' plot(ftrd.noip, model="Trended", type="hist")
#' plot(ftrd.noip, model="Trended", type="qq")
#' # dev.off()
#' 
nb_gof_m = function(counts, x, lib.sizes=colSums(counts), sim=999, model.fit="NBP", 
                     min.n=100, prior.df = 10, design = "simple"){
  
  stopifnot(model.fit %in% c("NBP","edgeR-common","edgeR-genewise","edgeR-trended",
                             "edgeR-tagcom","edgeR-tagtrd"))
  stopifnot(design %in% c("simple", "multiple"))
  
  # We recommend that genes with all zero counts be removed in advance;
  # The first argument, counts, should be already subsetted before passing to the function
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  N = m * n
  counts.dim = paste(m,"x",n)
  
  ## initialize simulation variables
  ord.res.sim.mat = matrix(0, nrow = (sim+1), ncol = N)   # ordered residual "big" matrix
  phi.hat.sim.mat = matrix(0, nrow = (sim+1), ncol = m)   # phi hat "big" matrix
  
  #### -----------------------------------------------------------------
  if (model.fit == "NBP"){
    mnbp.0 = model_nbp_m(counts, x, lib.sizes=colSums(counts))
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
#################################################################################
    ## the new method is applied here
 #  mnbp.h = model_nbp_m(y.mat.h, x, lib.sizes=colSums(y.mat.h))
      mnbp.h <- model_nbp_m_1(y.mat.h, x, lib.sizes=colSums(y.mat.h), obs.fit = mnbp.0)
#################################################################################
      ord.res.sim.mat[i, ] = mnbp.h$ord.res.vec
      phi.hat.sim.mat[i, ] = mnbp.h$phi.hat.mat[ ,1]
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
    phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0[ ,1]
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-common"){
    mcom.0 = model_edgeR_common(counts, x, lib.sizes=colSums(counts), design=design)
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
      mcom.h = model_edgeR_common(y.mat.h, x, lib.sizes=colSums(y.mat.h), design=design)
      ord.res.sim.mat[i, ] = mcom.h$ord.res.vec
      phi.hat.sim.mat[i, ] = mcom.h$phi.hat.mat
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
    phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-genewise"){
    mgen.0 = model_edgeR_genewise(counts, x, lib.sizes=colSums(counts), design=design)
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
      mgen.h = model_edgeR_genewise(y.mat.h, x, lib.sizes=colSums(y.mat.h), design=design)
      ord.res.sim.mat[i, ] = mgen.h$ord.res.vec
      phi.hat.sim.mat[i, ] = mgen.h$phi.hat.mat
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
    phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-trended"){
    mtrd.0 = model_edgeR_trended(counts, x, lib.sizes=colSums(counts), 
                                 min.n=min.n, design=design)
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
      mtrd.h = model_edgeR_trended(y.mat.h, x, lib.sizes=colSums(y.mat.h), 
                                   min.n=min.n, design=design)
      ord.res.sim.mat[i, ] = mtrd.h$ord.res.vec
      phi.hat.sim.mat[i, ] = mtrd.h$phi.hat.mat
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
    phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-tagcom"){
    mtgc.0 = model_edgeR_tagcom(counts, x, lib.sizes=colSums(counts), 
                                 prior.df = prior.df, design=design)
    mu.hat.mat0 = mtgc.0$mu.hat.mat
    phi.hat.mat0 = mtgc.0$phi.hat.mat
    res.omat0 = mtgc.0$res.omat
    ord.res.vec0 = mtgc.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      mtgc.h = model_edgeR_tagcom(y.mat.h, x, lib.sizes=colSums(y.mat.h), 
                                   prior.df = prior.df, design=design)
      ord.res.sim.mat[i, ] = mtgc.h$ord.res.vec
      phi.hat.sim.mat[i, ] = mtgc.h$phi.hat.mat
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
    phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  #### -----------------------------------------------------------------
  if (model.fit == "edgeR-tagtrd"){
    mtgt.0 = model_edgeR_tagtrd(counts, x, lib.sizes=colSums(counts), min.n=min.n,
                                 prior.df = prior.df, design=design)
    mu.hat.mat0 = mtgt.0$mu.hat.mat
    phi.hat.mat0 = mtgt.0$phi.hat.mat
    res.omat0 = mtgt.0$res.omat
    ord.res.vec0 = mtgt.0$ord.res.vec
    ## simulate new datasets and re-fit
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      mtgt.h = model_edgeR_tagtrd(y.mat.h, x, lib.sizes=colSums(y.mat.h), min.n=min.n,
                                   prior.df = prior.df, design=design)
      ord.res.sim.mat[i, ] = mtgt.h$ord.res.vec
      phi.hat.sim.mat[i, ] = mtgt.h$phi.hat.mat
    }
    close(pb)
    ord.res.sim.mat[(sim+1), ] = ord.res.vec0
    phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }

  #### -----------------------------------------------------------------
  # find the median of the big residual matrix (ordered): a vector
  ord.typ.res.sim = apply(ord.res.sim.mat[1:sim, ], 2, median)  # on simulated datasets ONLY!
  # subtract the typical residual vector from each row of the ordered residual matrix
  # dists.mat.res =  abs(sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))
  dists.mat.res = (sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))^2
  
  # construct new distance matrix D of dimension (R+1)-by-m
  grp.vec = ( seq_len( ncol(dists.mat.res) ) - 1 ) %/% n     # grouping vector
  dist.mat = t( rowsum(t(dists.mat.res), grp.vec) )    # vertical distance matrix (sim. + obs.)
  # THIS dist.mat IS UN-SORTED!! WE CAN USE THIS MATRIX FOR THE ENVELOPE METHOD CALCULATIONS
  pear.mat = t( rowsum(t(ord.res.sim.mat)^2, grp.vec) )  # Pearson stats matrix (sim. + obs.)
  
  #### -----------------------------------------------------------------
  ## calcualte test statistics and p-values for Monte Carlo method (sum of each row of dist.mat)
#   stat0.Vert = sum(dist.mat[(sim+1), ])  # squared vertical distance
#   stat0.Pear = sum(pear.mat[(sim+1), ])  # Pearson statistic
  
#   for (i in 1:sim){
#     stat.sim.Vert[i] = sum(dist.mat[i, ])
#     stat.sim.Pear[i] = sum(pear.mat[i, ])
#   }
#   pval.Vert = (sum(stat.sim.Vert >= stat0.Vert) + 1) / (sim + 1)      # ONE-SIDED!!
#   pval.Pear = (sum(stat.sim.Pear >= stat0.Pear) + 1) / (sim + 1)      # ONE-SIDED!!
#   pval.Pear = 2 * min( (sum(stat.sim.Pear >= stat0.Pear) + 1 ) / (sim + 1),
#                        1 - ( sum(stat.sim.Pear >= stat0.Pear) + 1 ) / (sim + 1) )  # TWO-SIDED
#   pv.Vert = round(pval.Vert, 6)  
#   pv.Pear = round(pval.Pear, 6)
  
  #### -----------------------------------------------------------------
  ## calcualte one p-value for each gene based on dist.mat or pear.mat, so a total of m p-values
  ## this p-value is the Monte Carlo p-value simply from a single univariate case
  v.pvals = numeric(m)     # M.C. p-values based on vertical distances
  p.pvals.1s = numeric(m)  # M.C. p-values based on Pearson statistics (1-sided)
  p.pvals.2s = numeric(m)  # M.C. p-values based on Pearson statistics (2-sided)
  for (i in 1:m){
    v.pvals[i] = (sum(dist.mat[1:sim,i] >= dist.mat[(sim+1),i]) + 1) / (sim + 1)
    p.pvals.1s[i] = (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1)
    p.pvals.2s[i] = 2 * min( (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1),
                          1 - (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1) )
  }
  
  # to avoid potential zero p-values
  v.pvals[v.pvals == 0] = 1/sim
  p.pvals.1s[p.pvals.1s == 0] = 1/sim
  p.pvals.2s[p.pvals.2s == 0] = 2/sim
  
  #### -----------------------------------------------------------------
  ## save as a list
  gof.obj = list(model.fit = model.fit,
                 counts.dim = counts.dim,
                 design.mat = x,
                 lib.sizes = lib.sizes,
                 v.pvals = v.pvals,
                 p.pvals.1s = p.pvals.1s,
                 p.pvals.2s = p.pvals.2s,
                 #mu.hat.m0 = mu.hat.mat0,
                 #ord.dist.mat = ord.dist.mat,
                 #ord.typ.dist  = ord.typ.dist,
                 #dist.obs = dist.obs,
                 phi.hat.sim.mat = phi.hat.sim.mat,
                 dist.mat = dist.mat,
                 pear.mat = pear.mat,
                 sim = sim
                 )
  
  # save the object as a "gofm" class
  class(gof.obj) = "gofm"
  return(gof.obj)
}



# #' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on an 
# #' RNA-Seq Dataset
# #' 
# #' @description This function performs simulation-based goodness-of-fit tests of different 
# #' negative binomial dispersion models for an RNA-Seq dataset with a known design 
# #' matrix. Currently supported models are implemented in the \code{NBPSeq} and \code{edgeR}
# #' packages.
# #' 
# #' @param counts an m-by-n count matrix of non-negative integers. For a typical
# #' RNA-Seq experiment, this is the read counts with m genes and n samples.
# #' @param x an n-by-p design matrix.
# #' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
# #' sums of the \code{counts} matrix.
# #' @param sim number of simulations performed.
# #' @param model.fit a string of characters specifying the negative binomial
# #' dispersion model used to fit the data. Currently the following dispersion models
# #' are available to be checked for goodness-of-fit:
# #' \itemize{
# #' \item NBP dispersion model in the \code{NBPSeq} package (\code{NBP})
# #' \item NB common dispersion model in the \code{\link{edgeR}} package (\code{edgeR-common})
# #' \item NB genewise dispersion model in the \code{\link{edgeR}} package (\code{edgeR-genewise})
# #' \item NB trended (non-parametric) dispersion model in the \code{\link{edgeR}} package
# #' (\code{edgeR-trended})
# #' \item NB tagwise-common dispersion model in the \code{\link{edgeR}} package 
# #' (\code{edgeR-tagcom})
# #' \item NB tagwise-trended dispersion model in the \code{\link{edgeR}} package 
# #' (\code{edgeR-tagtrd})
# #' }
# #' Users should specify \strong{exactly} the same characters as the
# #' ones in paratheses above for each dispersion model.
# #' @param min.n for \code{edgeR-trended} model only: specify the minimim number of genes in a bin
# #' (default is 10).
# #' @param prior.df control of the shrinkage applied for \code{edgeR} genewise, tagwise 
# #' (and its variants) dispersion models. Setting \code{prior.df=0} gives the genewise model
# #' without any shrinkage. Setting \code{prior.df=10} (default in \code{edgeR}) gives the tagwise
# #' model with shrinkage, either towards the common dispersion (\code{edgeR-tagcom}) or the
# #' global dispersion trend (\code{edgeR-tagtrd}).
# #' @param design specifications of testing (1) a single-group model (\code{single}); (2)
# #' a multiple-group model (\code{multiple}); (3) complex design with interactions
# #' (\code{complex}). \strong{The codes are only tested for single- and multiple-group cases, and
# #' the modeling of dispersions in the multiple-group case is still under consideration!
# #' For the complex design with interactions, though we can pass the design matrix directly
# #' in edgeR, the results may not make any sense!}
# #' 
# #' @return An object of class "gofm" to which other methods can be applied.
# #' 
# #' @details When the response is a count matrix, we can use this function to test
# #' the goodness-of-fit of a specified negative binomial dispersion model. Specifically,
# #' this function calls \code{\link{model_nbp_m}} to fit the NBP dispersion model, calls
# #' \code{\link{model_edgeR_common}} to fit the common dispersion model, calls
# #' \code{\link{model_edgeR_genewise}} to fit the pure genewise dispersion model (without
# #' shrinkage), calls \code{\link{model_edgeR_trended}} to fit the non-parametric trended
# #' dispersion model in \code{edgeR}, calls \code{\link{model_edgeR_tagcom}} to fit the tagwise
# #' model with the dispersions shrinking towards a common dispersion, and calls
# #' \code{\link{model_edgeR_tagtrd}} to fit the tagwise model with the dispersions shrinking
# #' towards the global dispersion trend.
# #' 
# #' @usage 
# #' nb_gof_m(counts, x, lib.sizes=colSums(counts), sim=999, model.fit="NBP", min.n=100, 
# #' prior.df = 10, design = "simple")
# #' 
# #' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
# #' 
# #' @export
# #' 
# #' @references 
# #' Mi G, Di Y and Schafer DW (2013). Goodness-of-Fit Tests and Model Diagnostics
# #' for Negative Binomial Regression of RNA sequencing Data. \emph{Biometrics} (submitted).
# #' 
# #' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
# #' Model for Assessing Differential Gene Expression from RNA-Seq", \emph{Statistical 
# #' Applications in Genetics and Molecular Biology}, 10 (1).
# #' 
# #' McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression analysis of multifactor 
# #' RNA-Seq experiments with respect to biological variation. \emph{Nucleic Acids Research} 
# #' 40, 4288-4297.
# #' 
# #' Cox, DR, and Reid, N (1987). Parameter orthogonality and approximate conditional inference.
# #' \emph{Journal of the Royal Statistical Society Series B} 49, 1-39.
# #' 
# #' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
# #' 
# #' @examples 
# #' ## load package into session:
# #' library(NBGOF)
# #' 
# #' ## basic set-up of the model:
# #' seed = 315825
# #' sim = 999
# #' conf.env = 0.95
# #' 
# #' ## basic parameter specifications:
# #' m = 1000; # 1000 genes
# #' n = 6;    # 6 samples
# #' s = 1e6;  # lib.sizes
# #' offset = log(s);  # make sure offset for NBP functions should take the log
# #' offset    # 13.82
# #' 
# #' ## simulate coefficients:
# #' beta = matrix(0, m, 1);  # beta0 only
# #' set.seed(seed)
# #' beta[,1] = rnorm(m, 5, 2) - offset;   # beta0 ~ normal (mean and a based on real RNA-Seq data)
# #' beta  # m-by-1 matrix
# #' 
# #' ## design matrix (single group only):
# #' x = matrix(rep(1,n))
# #' x
# #' 
# #' ## specify mean levels:
# #' mu = round(t(s * exp(x %*% t(beta))));
# #' sum(rowSums(mu) == 0)   # 4
# #' mu[rowSums(mu) == 0, ] = 1
# #' sum(rowSums(mu) == 0)   # 0
# #' pi = mu/s
# #' 
# #' ## simulate an m-by-n count matrix mimicking a RNA-Seq dataset:
# #' ## simulate NB1.8 data with random uniform noise added to the dispersion
# #' alpha1 = -0.2  # NB1.8
# #' phi0 = 0.02
# #' alpha0 = log(phi0)
# #' phi.nbp = phi0 * pi^alpha1
# #' range(phi.nbp)
# #' cbind(mu[,1], phi.nbp[,1]);
# #' 
# #' ## add noise:
# #' a = 0.1
# #' set.seed(seed)  
# #' phi.noi.vec = phi.nbp[ ,1] * exp(matrix(runif(m, -a, a), nr=m, nc=1))
# #' phi.noi = matrix(phi.noi.vec, nr=m, nc=n)
# #' range(phi.noi)
# #' cbind(mu[,1], phi.noi[,1])  # make sure phi's are in reasonable range
# #' 
# #' ## generate NBP response with added noise:
# #' set.seed(seed)
# #' y = rnbinom(m * n, mu=mu, size=1/phi.noi)
# #' dim(y) = dim(mu)
# #' rownames(y) = paste("g", seq(1,m), sep="")
# #' colnames(y) = paste("s", seq(1,n), sep="")
# #' plot(mu, phi.noi, log="xy")
# #' 
# #' ## GOF tests for different dispersion models:
# #' ## CAUTION: may be time-consuming depending on the size of data and simulations
# #' fnbp.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="NBP")
# #' fcom.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-common")
# #' fgen.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-genewise")
# #' ftgc.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-tagcom")
# #' ftgt.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-tagtrd")
# #' ftrd.noip = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-trended")
# #' 
# #' ## summarize the GOF test results:
# #' summary(fnbp.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
# #' summary(fcom.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
# #' summary(fgen.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
# #' summary(ftgc.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
# #' summary(ftgt.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
# #' summary(ftrd.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
# #' 
# #' ## evaluate GOF p-values using histograms and quantile-quantile uniform plots:
# #' # pdf(file="plots-NBP-01-1grp.pdf", width=8, height=24)
# #' par(mfrow=c(6,2))
# #' plot(fnbp.noip, model="NBP", type="hist")
# #' plot(fnbp.noip, model="NBP", type="qq")
# #' plot(fcom.noip, model="Common", type="hist")
# #' plot(fcom.noip, model="Common", type="qq")
# #' plot(fgen.noip, model="Genewise", type="hist")
# #' plot(fgen.noip, model="Genewise", type="qq")
# #' plot(ftgc.noip, model="Tagwise-Common", type="hist")
# #' plot(ftgc.noip, model="Tagwise-Common", type="qq")
# #' plot(ftgt.noip, model="Tagwise-Trended", type="hist")
# #' plot(ftgt.noip, model="Tagwise-Trended", type="qq")
# #' plot(ftrd.noip, model="Trended", type="hist")
# #' plot(ftrd.noip, model="Trended", type="qq")
# #' # dev.off()
# #' 
# nb_gof_m = function(counts, x, lib.sizes=colSums(counts), sim=999, model.fit="NBP", 
#                      min.n=100, prior.df = 10, design = "simple"){
#   
#   stopifnot(model.fit %in% c("NBP","edgeR-common","edgeR-genewise","edgeR-trended",
#                              "edgeR-tagcom","edgeR-tagtrd"))
#   stopifnot(design %in% c("simple", "multiple"))
#   
#   # We recommend that genes with all zero counts be removed in advance;
#   # The first argument, counts, should be already subsetted before passing to the function
#   
#   m = dim(counts)[1]
#   n = dim(counts)[2]
#   N = m * n
#   counts.dim = paste(m,"x",n)
#   
#   ## initialize simulation variables
#   ord.res.sim.mat = matrix(0, nrow = (sim+1), ncol = N)   # ordered residual "big" matrix
#   phi.hat.sim.mat = matrix(0, nrow = (sim+1), ncol = m)   # phi hat "big" matrix
#   
#   #### -----------------------------------------------------------------
#   if (model.fit == "NBP"){
#     mnbp.0 = model_nbp_m(counts, x, lib.sizes=colSums(counts))
#     mu.hat.mat0 = mnbp.0$mu.hat.mat
#     phi.hat.mat0 = mnbp.0$phi.hat.mat
#     res.omat0 = mnbp.0$res.omat
#     ord.res.vec0 = mnbp.0$ord.res.vec
#     ## simulate new datasets and re-fit
#     pb = txtProgressBar(style=3)
#     ## the new method is applied here
#     for (i in 1:sim){
#       setTxtProgressBar(pb, i/sim)
#       y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
#       dim(y.mat.h) = dim(counts)
#       rownames(y.mat.h) = rownames(counts)
#       colnames(y.mat.h) = colnames(counts)
#       mnbp.h = model_nbp_m_1(y.mat.h, x, lib.sizes=colSums(y.mat.h))
#       ord.res.sim.mat[i, ] = mnbp.h$ord.res.vec
#       phi.hat.sim.mat[i, ] = mnbp.h$phi.hat.mat[ ,1]
#     }
#     close(pb)
#     ord.res.sim.mat[(sim+1), ] = ord.res.vec0
#     phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0[ ,1]
#   }
#   #### -----------------------------------------------------------------
#   
#   if (model.fit == "edgeR-common"){
#     mcom.0 = model_edgeR_common(counts, x, lib.sizes=colSums(counts), design=design)
#     mu.hat.mat0 = mcom.0$mu.hat.mat
#     phi.hat.mat0 = mcom.0$phi.hat.mat
#     res.omat0 = mcom.0$res.omat
#     ord.res.vec0 = mcom.0$ord.res.vec
#     ## simulate new datasets and re-fit
#     pb = txtProgressBar(style=3)
#     for (i in 1:sim){
#       setTxtProgressBar(pb, i/sim)
#       y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
#       dim(y.mat.h) = dim(counts)
#       rownames(y.mat.h) = rownames(counts)
#       colnames(y.mat.h) = colnames(counts)
#       mcom.h = model_edgeR_common(y.mat.h, x, lib.sizes=colSums(y.mat.h), design=design)
#       ord.res.sim.mat[i, ] = mcom.h$ord.res.vec
#       phi.hat.sim.mat[i, ] = mcom.h$phi.hat.mat
#     }
#     close(pb)
#     ord.res.sim.mat[(sim+1), ] = ord.res.vec0
#     phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
#   }
#   #### -----------------------------------------------------------------
#   if (model.fit == "edgeR-genewise"){
#     mgen.0 = model_edgeR_genewise(counts, x, lib.sizes=colSums(counts), design=design)
#     mu.hat.mat0 = mgen.0$mu.hat.mat
#     phi.hat.mat0 = mgen.0$phi.hat.mat
#     res.omat0 = mgen.0$res.omat
#     ord.res.vec0 = mgen.0$ord.res.vec
#     ## simulate new datasets and re-fit
#     pb = txtProgressBar(style=3)
#     for (i in 1:sim){
#       setTxtProgressBar(pb, i/sim)
#       y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
#       dim(y.mat.h) = dim(counts)
#       rownames(y.mat.h) = rownames(counts)
#       colnames(y.mat.h) = colnames(counts)
#       mgen.h = model_edgeR_genewise(y.mat.h, x, lib.sizes=colSums(y.mat.h), design=design)
#       ord.res.sim.mat[i, ] = mgen.h$ord.res.vec
#       phi.hat.sim.mat[i, ] = mgen.h$phi.hat.mat
#     }
#     close(pb)
#     ord.res.sim.mat[(sim+1), ] = ord.res.vec0
#     phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
#   }
#   #### -----------------------------------------------------------------
#   if (model.fit == "edgeR-trended"){
#     mtrd.0 = model_edgeR_trended(counts, x, lib.sizes=colSums(counts), 
#                                  min.n=min.n, design=design)
#     mu.hat.mat0 = mtrd.0$mu.hat.mat
#     phi.hat.mat0 = mtrd.0$phi.hat.mat
#     res.omat0 = mtrd.0$res.omat
#     ord.res.vec0 = mtrd.0$ord.res.vec
#     ## simulate new datasets and re-fit
#     pb = txtProgressBar(style=3)
#     for (i in 1:sim){
#       setTxtProgressBar(pb, i/sim)
#       y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
#       dim(y.mat.h) = dim(counts)
#       rownames(y.mat.h) = rownames(counts)
#       colnames(y.mat.h) = colnames(counts)
#       mtrd.h = model_edgeR_trended(y.mat.h, x, lib.sizes=colSums(y.mat.h), 
#                                    min.n=min.n, design=design)
#       ord.res.sim.mat[i, ] = mtrd.h$ord.res.vec
#       phi.hat.sim.mat[i, ] = mtrd.h$phi.hat.mat
#     }
#     close(pb)
#     ord.res.sim.mat[(sim+1), ] = ord.res.vec0
#     phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
#   }
#   #### -----------------------------------------------------------------
#   if (model.fit == "edgeR-tagcom"){
#     mtgc.0 = model_edgeR_tagcom(counts, x, lib.sizes=colSums(counts), 
#                                  prior.df = prior.df, design=design)
#     mu.hat.mat0 = mtgc.0$mu.hat.mat
#     phi.hat.mat0 = mtgc.0$phi.hat.mat
#     res.omat0 = mtgc.0$res.omat
#     ord.res.vec0 = mtgc.0$ord.res.vec
#     ## simulate new datasets and re-fit
#     pb = txtProgressBar(style=3)
#     for (i in 1:sim){
#       setTxtProgressBar(pb, i/sim)
#       y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
#       dim(y.mat.h) = dim(counts)
#       rownames(y.mat.h) = rownames(counts)
#       colnames(y.mat.h) = colnames(counts)
#       mtgc.h = model_edgeR_tagcom(y.mat.h, x, lib.sizes=colSums(y.mat.h), 
#                                    prior.df = prior.df, design=design)
#       ord.res.sim.mat[i, ] = mtgc.h$ord.res.vec
#       phi.hat.sim.mat[i, ] = mtgc.h$phi.hat.mat
#     }
#     close(pb)
#     ord.res.sim.mat[(sim+1), ] = ord.res.vec0
#     phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
#   }
#   #### -----------------------------------------------------------------
#   if (model.fit == "edgeR-tagtrd"){
#     mtgt.0 = model_edgeR_tagtrd(counts, x, lib.sizes=colSums(counts), min.n=min.n,
#                                  prior.df = prior.df, design=design)
#     mu.hat.mat0 = mtgt.0$mu.hat.mat
#     phi.hat.mat0 = mtgt.0$phi.hat.mat
#     res.omat0 = mtgt.0$res.omat
#     ord.res.vec0 = mtgt.0$ord.res.vec
#     ## simulate new datasets and re-fit
#     pb = txtProgressBar(style=3)
#     for (i in 1:sim){
#       setTxtProgressBar(pb, i/sim)
#       y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
#       dim(y.mat.h) = dim(counts)
#       rownames(y.mat.h) = rownames(counts)
#       colnames(y.mat.h) = colnames(counts)
#       mtgt.h = model_edgeR_tagtrd(y.mat.h, x, lib.sizes=colSums(y.mat.h), min.n=min.n,
#                                    prior.df = prior.df, design=design)
#       ord.res.sim.mat[i, ] = mtgt.h$ord.res.vec
#       phi.hat.sim.mat[i, ] = mtgt.h$phi.hat.mat
#     }
#     close(pb)
#     ord.res.sim.mat[(sim+1), ] = ord.res.vec0
#     phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
#   }
# 
#   #### -----------------------------------------------------------------
#   # find the median of the big residual matrix (ordered): a vector
#   ord.typ.res.sim = apply(ord.res.sim.mat[1:sim, ], 2, median)  # on simulated datasets ONLY!
#   # subtract the typical residual vector from each row of the ordered residual matrix
#   # dists.mat.res =  abs(sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))
#   dists.mat.res = (sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))^2
#   
#   # construct new distance matrix D of dimension (R+1)-by-m
#   grp.vec = ( seq_len( ncol(dists.mat.res) ) - 1 ) %/% n     # grouping vector
#   dist.mat = t( rowsum(t(dists.mat.res), grp.vec) )    # vertical distance matrix (sim. + obs.)
#   # THIS dist.mat IS UN-SORTED!! WE CAN USE THIS MATRIX FOR THE ENVELOPE METHOD CALCULATIONS
#   pear.mat = t( rowsum(t(ord.res.sim.mat)^2, grp.vec) )  # Pearson stats matrix (sim. + obs.)
#   
#   #### -----------------------------------------------------------------
#   ## calcualte test statistics and p-values for Monte Carlo method (sum of each row of dist.mat)
# #   stat0.Vert = sum(dist.mat[(sim+1), ])  # squared vertical distance
# #   stat0.Pear = sum(pear.mat[(sim+1), ])  # Pearson statistic
#   
# #   for (i in 1:sim){
# #     stat.sim.Vert[i] = sum(dist.mat[i, ])
# #     stat.sim.Pear[i] = sum(pear.mat[i, ])
# #   }
# #   pval.Vert = (sum(stat.sim.Vert >= stat0.Vert) + 1) / (sim + 1)      # ONE-SIDED!!
# #   pval.Pear = (sum(stat.sim.Pear >= stat0.Pear) + 1) / (sim + 1)      # ONE-SIDED!!
# #   pval.Pear = 2 * min( (sum(stat.sim.Pear >= stat0.Pear) + 1 ) / (sim + 1),
# #                        1 - ( sum(stat.sim.Pear >= stat0.Pear) + 1 ) / (sim + 1) )  # TWO-SIDED
# #   pv.Vert = round(pval.Vert, 6)  
# #   pv.Pear = round(pval.Pear, 6)
#   
#   #### -----------------------------------------------------------------
#   ## calcualte one p-value for each gene based on dist.mat or pear.mat, so a total of m p-values
#   ## this p-value is the Monte Carlo p-value simply from a single univariate case
#   v.pvals = numeric(m)     # M.C. p-values based on vertical distances
#   p.pvals.1s = numeric(m)  # M.C. p-values based on Pearson statistics (1-sided)
#   p.pvals.2s = numeric(m)  # M.C. p-values based on Pearson statistics (2-sided)
#   for (i in 1:m){
#     v.pvals[i] = (sum(dist.mat[1:sim,i] >= dist.mat[(sim+1),i]) + 1) / (sim + 1)
#     p.pvals.1s[i] = (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1)
#     p.pvals.2s[i] = 2 * min( (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1),
#                           1 - (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1) )
#   }
#   
#   # to avoid potential zero p-values
#   v.pvals[v.pvals == 0] = 1/sim
#   p.pvals.1s[p.pvals.1s == 0] = 1/sim
#   p.pvals.2s[p.pvals.2s == 0] = 2/sim
#   
#   #### -----------------------------------------------------------------
#   ## save as a list
#   gof.obj = list(model.fit = model.fit,
#                  counts.dim = counts.dim,
#                  design.mat = x,
#                  lib.sizes = lib.sizes,
#                  v.pvals = v.pvals,
#                  p.pvals.1s = p.pvals.1s,
#                  p.pvals.2s = p.pvals.2s,
#                  #mu.hat.m0 = mu.hat.mat0,
#                  #ord.dist.mat = ord.dist.mat,
#                  #ord.typ.dist  = ord.typ.dist,
#                  #dist.obs = dist.obs,
#                  phi.hat.sim.mat = phi.hat.sim.mat,
#                  dist.mat = dist.mat,
#                  pear.mat = pear.mat,
#                  sim = sim
#                  )
#   
#   # save the object as a "gofm" class
#   class(gof.obj) = "gofm"
#   return(gof.obj)
# }
