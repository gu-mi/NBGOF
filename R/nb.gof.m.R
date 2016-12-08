
#' @title Main Function of Implementing Simulation-based Goodness-of-Fit Tests on an RNA-Seq Dataset
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

#' @param model a string of characters specifying the negative binomial
#' dispersion model used to fit the data. Currently the following dispersion models
#' are available to be checked for goodness-of-fit:
#' \itemize{
#' \item NBP dispersion model in the \code{NBPSeq} package (\code{NBP})
#' \item NBQ dispersion model in the \code{NBPSeq} pacakge (\code{NBQ})
#' \item NBS dispersion model in the \code{NBPSeq} pacakge (\code{NBS})
#' \item NBSTEP dispersion model in the \code{NBPSeq} pacakge (\code{NBSTEP})
#' \item NB common dispersion model in the \code{edgeR} package (\code{Common})
#' \item NB genewise in the \code{NBPSeq} package (\code{Genewise})
#' \item NB trended (non-parametric) dispersion model in the \code{edgeR} package
#' (\code{Trended})
#' \item NB tagwise-common dispersion model in the \code{edgeR} package 
#' (\code{Tagwise-Common})
#' \item NB tagwise-trended dispersion model in the \code{edgeR} package 
#' (\code{Tagwise-Trend})
#' }
#' Users should specify \strong{exactly} the same characters as the
#' ones in paratheses above for each dispersion model.

#' @param method method for estimating dispersions. MAPL: maximum adjusted profile likelihood other estimation methods from
#' \code{edgeR} can also be specified.

#' @param min.n for \code{Trended} model only: specify the minimim number of genes in a bin (default is 100).

#' @param prior.df control of the shrinkage applied for \code{edgeR} genewise, tagwise 
#' (and its variants) dispersion models. Setting \code{prior.df=0} gives the genewise model
#' without any shrinkage. Setting \code{prior.df=10} (default in \code{edgeR}) gives the tagwise
#' model with shrinkage, either towards the common dispersion (\code{Tagwise-Common}) or the
#' global dispersion trend (\code{Tagwise-Trend}). In the paper of McCarthy, DJ et al. (2012),
#' in the formula G_0 = prior.df/df = 20/df: df is the "residual df" (equals the number of libraries minus the number of distinct
#' treatment groups), and prior.df=20. For example, if there are six libraries and two distinct groups, 
#' df = residual df = 4. If we set prior.df = 20, then G_0 = 5. G_0 is the "prior.n" in \code{edgeR}'s earlier versions.
#' Users should be careful in choosing appropriate prior.df when using tagwise models.

#' @param seed initial random number seed for reproducibility of re-simulations (default is 1).

#' @param ncores number of CPU cores to use. If unspecified, use the total number of CPUs minus 1.
#' 
#' @param ... (for future use)
#' 
#' @return An object of class "gofm" to which other methods can be applied.
#' 
#' @details When the response is a count matrix, we can use this function to test
#' the goodness-of-fit of a specified negative binomial dispersion model. Specifically,
#' this function calls \code{\link{model.nbp.m}} to fit the NBP dispersion model, 
#' calls \code{\link{model.nbq.m}} to fit the NBQ dispersion model, calls
#' \code{\link{model.edgeR.common}} to fit the common dispersion model, calls
#' \code{\link{model.edgeR.genewise}} to fit the pure genewise dispersion model (without
#' shrinkage), calls \code{\link{model.edgeR.trended}} to fit the non-parametric trended
#' dispersion model in \code{edgeR}, calls \code{\link{model.edgeR.tagcom}} to fit the tagwise
#' model with the dispersions shrinking towards a common dispersion, and calls
#' \code{\link{model.edgeR.tagtrd}} to fit the tagwise model with the dispersions shrinking
#' towards the global dispersion trend.
#' 
#' @usage 
#' nb.gof.m(counts, x, lib.sizes=colSums(counts), sim=999, model=NULL, method=NULL, 
#' min.n=100, prior.df = 10, seed=1, ncores = NULL, ...)
#' 
#' @author Gu Mi <neo.migu@gmail.com>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references 
#' Mi, G, Di, Y, & Schafer, DW (2015). Goodness-of-Fit Tests and Model Diagnostics for 
#' Negative Binomial Regression of RNA Sequencing Data. \emph{PLOS ONE}, 10(3).
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
#' seed = 539
#' sim = 999
#' conf.env = 0.95
#' 
#' ## basic parameter specifications:
#' m = 1000 # 1000 genes
#' n = 6    # 6 samples
#' s = 1e6  # lib.sizes
#' offset = log(s)  # make sure offset for NBP functions should take the log
#' 
#' ## simulate coefficients:
#' beta = matrix(0, m, 1)  # beta0 only
#' set.seed(seed)
#' beta[,1] = rnorm(m, 5, 2) - offset   # beta0 ~ normal (mean and a based on real RNA-Seq data)
#' 
#' ## design matrix (single group only):
#' x = matrix(rep(1,n))
#' 
#' ## specify mean levels:
#' mu = round(t(s * exp(x %*% t(beta))))
#' pi = mu/s
#' 
#' ## simulate an m-by-n count matrix mimicking a RNA-Seq dataset:
#' ## simulate NB1.8 data with random uniform noise added to the dispersion
#' alpha1 = -0.2  # NB1.8
#' phi0 = 0.02
#' alpha0 = log(phi0)
#' phi.nbp = phi0 * pi^alpha1
#' range(phi.nbp)
#' cbind(mu[,1], phi.nbp[,1])
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
#' ## GOF tests for different dispersion models, using parallel computing:
#' ## CAUTION: may be time-consuming depending on the size of data and simulations!
#' nc = detectCores() - 1
#' fnbp.noip = nb.gof.m(counts=y, x=x, sim=sim, model="NBP", method="MAPL", seed=1, ncores=nc)
#' fnbq.noip = nb.gof.m(counts=y, x=x, sim=sim, model="NBQ", method="MAPL", seed=1, ncores=nc)
#' fcom.noip = nb.gof.m(counts=y, x=x, sim=sim, model="Common", method="CoxReid", seed=1, ncores=nc)
#' fgen.noip = nb.gof.m(counts=y, x=x, sim=sim, model="Genewise", method="auto", seed=1, ncores=nc)
#' ftgc.noip = nb.gof.m(counts=y, x=x, sim=sim, model="Tagwise-Common", method="CoxReid", seed=1, ncores=nc)
#' ftgt.noip = nb.gof.m(counts=y, x=x, sim=sim, model="Tagwise-Trend", method="auto", seed=1, ncores=nc)
#' ftrd.noip = nb.gof.m(counts=y, x=x, sim=sim, model="Trended", method="auto", seed=1, ncores=nc)
#' 
#' ## summarize the GOF test results:
#' summary(fnbp.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(fnbq.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(fcom.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(fgen.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(ftgc.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(ftgt.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' summary(ftrd.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
#' 
#' ## evaluate GOF p-values using histograms and quantile-quantile uniform plots:
#' # pdf("plots-NBP-01-1grp.pdf", width=20, height=20)
#' h1 = diagPlot(fnbp.noip, type="hist", binwidth=0.1)
#' q1 = diagPlot(fnbp.noip, type="qq")
#' h2 = diagPlot(fnbq.noip, type="hist", binwidth=0.1)
#' q2 = diagPlot(fnbq.noip, type="qq")
#' h3 = diagPlot(fcom.noip, type="hist", binwidth=0.1)
#' q3 = diagPlot(fcom.noip, type="qq")
#' h4 = diagPlot(fgen.noip, type="hist", binwidth=0.1)
#' q4 = diagPlot(fgen.noip, type="qq")
#' h5 = diagPlot(ftgc.noip, type="hist", binwidth=0.1)
#' q5 = diagPlot(ftgc.noip, type="qq")
#' h6 = diagPlot(ftgt.noip, type="hist", binwidth=0.1)
#' q6 = diagPlot(ftgt.noip, type="qq")
#' h7 = diagPlot(ftrd.noip, type="hist", binwidth=0.1)
#' q7 = diagPlot(ftrd.noip, type="qq")
#' multiplot(h1, q1, h2, q2, h3, q3, h4, q4, h5, q5, h6, q6, h7, q7, cols=4, layout=matrix(seq(1,16), nrow=4, byrow=TRUE))
#' # dev.off()
#' 
nb.gof.m = function(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), sim=999, model=NULL, method = NULL, min.n=100, prior.df = 10, 
                    seed=1, ncores = NULL, ...){
  
#   # model specifications
#   if (is.null(model)){
#     stop("You must specify a model name!")
#   }
#   if ( !model %in% c("NBP", "NBQ", "NBS", "STEP", "Common", "Genewise", "Trended", "Tagwise-Common", "Tagwise-Trend") ) {
#     stop("You must specify one of these model names: 
#          'NBP', 'NBQ', 'NBS', 'STEP', 'Common', 'Genewise', 'Trended', 'Tagwise-Common' or 'Tagwise-Trend'!")
#   }
  
#   # method specifications
#   if ( model %in% c("NBP", "NBQ", "NBS", "STEP") & !method %in% c("ML", "MAPL") ) {
#     stop("You must specify one of these estimation methods (in NBPSeq): 'ML' or 'MAPL'!")
#   }
#   if ( model %in% c("Common", "Tagwise-Common") & !method %in% c("CoxReid", "Pearson", "deviance") ) {
#     stop("You must specify one of these methods for estimating the dispersion (in edgeR): 'CoxReid', 'Pearson' or 'deviance'!")
#   }
#   if ( model %in% c("Trended", "Tagwise-Trend", "Genewise") & !method %in% c("auto", "bin.spline", "bin.loess", "power", "spline") ) {
#     stop("You must specify one of these methods for estimating the trended dispersions (in edgeR): 'auto', 'bin.spline','bin.loess', 'power', or 'spline'!")
#   }
  
  # parallel computing: specify number of cores to use
  if (is.null(ncores)){
    ncores = detectCores() - 1
  }
  # register for foreach
  registerDoMC(ncores)
  
  # keep only complete cases:
  counts = counts[complete.cases(counts), ]
  
  # We recommend that genes with all zero counts be removed in advance
  # The first argument, counts, should be already subsetted before passing to the function
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  N = m * n
  counts.dim = paste(m,"x",n)
  
  ## initialize simulation variables
  #ord.res.sim.mat = matrix(0, nrow = (sim+1), ncol = N)   # ordered residual "big" matrix
  #phi.hat.sim.mat = matrix(0, nrow = (sim+1), ncol = m)   # phi hat "big" matrix
  
  #### -----------------------------------------------------------------
  if (model == "NBP"){
    mnbp.0 = model.nbp.m(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), method=method)
    mu.hat.mat0 = mnbp.0$mu.hat.mat
    phi.hat.mat0 = mnbp.0$phi.hat.mat
    res.omat0 = mnbp.0$res.omat
    ord.res.vec0 = mnbp.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.nbp.m(y.mat.h, x, lib.sizes=colSums(y.mat.h), method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0[ ,1]
  }
  
  #### -----------------------------------------------------------------
  if (model == "NBQ"){
    mnbq.0 = model.nbq.m(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), method=method)
    mu.hat.mat0 = mnbq.0$mu.hat.mat
    phi.hat.mat0 = mnbq.0$phi.hat.mat
    res.omat0 = mnbq.0$res.omat
    ord.res.vec0 = mnbq.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.nbq.m(y.mat.h, x, lib.sizes=colSums(y.mat.h), method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0[ ,1]
  }
  
  #### -----------------------------------------------------------------
  if (model == "NBS"){
    mnbs.0 = model.nbs.m(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), method=method)
    mu.hat.mat0 = mnbs.0$mu.hat.mat
    phi.hat.mat0 = mnbs.0$phi.hat.mat
    res.omat0 = mnbs.0$res.omat
    ord.res.vec0 = mnbs.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.nbs.m(y.mat.h, x, lib.sizes=colSums(y.mat.h), method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0[ ,1]
  }
  
  #### -----------------------------------------------------------------
  if (model == "STEP"){
    mnbstep.0 = model.nbstep.m(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), method=method)
    mu.hat.mat0 = mnbstep.0$mu.hat.mat
    phi.hat.mat0 = mnbstep.0$phi.hat.mat
    res.omat0 = mnbstep.0$res.omat
    ord.res.vec0 = mnbstep.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.nbstep.m(y.mat.h, x, lib.sizes=colSums(y.mat.h), method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0[ ,1]
  }
  
  #### -----------------------------------------------------------------
  if (model == "Common"){
    mcom.0 = model.edgeR.common(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), method=method)
    mu.hat.mat0 = mcom.0$mu.hat.mat
    phi.hat.mat0 = mcom.0$phi.hat.mat
    res.omat0 = mcom.0$res.omat
    ord.res.vec0 = mcom.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.edgeR.common(y.mat.h, x, lib.sizes=colSums(y.mat.h), method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  
  #### -----------------------------------------------------------------
  if (model == "Genewise"){
    mgen.0 = model.edgeR.genewise(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), min.n=min.n, method=method)
    mu.hat.mat0 = mgen.0$mu.hat.mat
    phi.hat.mat0 = mgen.0$phi.hat.mat
    res.omat0 = mgen.0$res.omat
    ord.res.vec0 = mgen.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.edgeR.genewise(y.mat.h, x, lib.sizes=colSums(y.mat.h), min.n=min.n, method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  
#   #### -----------------------------------------------------------------
    # genewise method in NBPSeq: TBD
  if (model == "Genewise-HOA"){
    mgen.0 = model.genewise(counts, x)
    mu.hat.mat0 = mgen.0$mu.hat.mat
    phi.hat.mat0 = mgen.0$phi.hat.mat
    res.omat0 = mgen.0$res.omat
    ord.res.vec0 = mgen.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.genewise(y.mat.h, x)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  
  #### -----------------------------------------------------------------
  if (model == "Trended"){
    mtrd.0 = model.edgeR.trended(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), min.n=min.n, method=method)
    mu.hat.mat0 = mtrd.0$mu.hat.mat
    phi.hat.mat0 = mtrd.0$phi.hat.mat
    res.omat0 = mtrd.0$res.omat
    ord.res.vec0 = mtrd.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.edgeR.trended(y.mat.h, x, lib.sizes=colSums(y.mat.h), min.n=min.n, method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  
  #### -----------------------------------------------------------------
  if (model == "Tagwise-Common"){
    mtgc.0 = model.edgeR.tagcom(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), prior.df = prior.df, method=method)
    mu.hat.mat0 = mtgc.0$mu.hat.mat
    phi.hat.mat0 = mtgc.0$phi.hat.mat
    res.omat0 = mtgc.0$res.omat
    ord.res.vec0 = mtgc.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.edgeR.tagcom(y.mat.h, x, lib.sizes=colSums(y.mat.h), prior.df = prior.df, method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }
  
  #### -----------------------------------------------------------------
  if (model == "Tagwise-Trend"){
    mtgt.0 = model.edgeR.tagtrd(counts, x, lib.sizes=colSums(counts, na.rm = TRUE), min.n=min.n, prior.df = prior.df, method=method)
    mu.hat.mat0 = mtgt.0$mu.hat.mat
    phi.hat.mat0 = mtgt.0$phi.hat.mat
    res.omat0 = mtgt.0$res.omat
    ord.res.vec0 = mtgt.0$ord.res.vec
    ## simulate new datasets and re-fit
    #pb = txtProgressBar(style=3)
    ord.res.sim.mat.tmp = foreach(i=1:sim, .combine="rbind", .inorder=TRUE) %dopar% {
      #setTxtProgressBar(pb, i/sim)
      set.seed(i+seed)
      y.mat.h = rnbinom(n=N, mu=mu.hat.mat0, size=1/phi.hat.mat0)
      dim(y.mat.h) = dim(counts)
      rownames(y.mat.h) = rownames(counts)
      colnames(y.mat.h) = colnames(counts)
      model.edgeR.tagtrd(y.mat.h, x, lib.sizes=colSums(y.mat.h), min.n=min.n, prior.df = prior.df, method=method)$ord.res.vec
    }
    #close(pb)
    dimnames(ord.res.sim.mat.tmp) = NULL
    ord.res.sim.mat = rbind(ord.res.sim.mat.tmp, ord.res.vec0)
    #phi.hat.sim.mat[(sim+1), ] = phi.hat.mat0
  }

  #### -----------------------------------------------------------------
  # find the median of the big residual matrix (ordered): a vector
  ord.typ.res.sim = apply(ord.res.sim.mat[1:sim, ], 2, median)  # on simulated datasets ONLY!
  # subtract the typical residual vector from each row of the ordered residual matrix
  dists.mat.res = (sweep(ord.res.sim.mat, 2, ord.typ.res.sim, "-"))^2
  
  # construct new distance matrix D of dimension (R+1)-by-m
  grp.vec = ( seq_len( ncol(dists.mat.res) ) - 1 ) %/% n     # grouping vector
  dist.mat = t( rowsum(t(dists.mat.res), grp.vec) )    # vertical distance matrix (sim. + obs.)
  # THIS dist.mat IS UN-SORTED
  pear.mat = t( rowsum(t(ord.res.sim.mat)^2, grp.vec) )  # Pearson stats matrix (sim. + obs.)
  
  #### -----------------------------------------------------------------
  ## calcualte one p-value for each gene based on dist.mat or pear.mat, so a total of m p-values
  ## this p-value is the Monte Carlo p-value simply from a single univariate case
  v.pvals = numeric(m)     # M.C. p-values based on vertical distances
  p.pvals.1s = numeric(m)  # M.C. p-values based on Pearson statistics (1-sided)
  for (i in 1:m){
    v.pvals[i] = (sum(dist.mat[1:sim,i] >= dist.mat[(sim+1),i]) + 1) / (sim + 1)
    p.pvals.1s[i] = (sum(pear.mat[1:sim,i] >= pear.mat[(sim+1),i]) + 1) / (sim + 1)
  }
  
  # to avoid potential zero p-values
  v.pvals[v.pvals == 0] = 1/sim
  p.pvals.1s[p.pvals.1s == 0] = 1/sim
  
  #### -----------------------------------------------------------------
  ## save as a list
  gof.obj = list(model = model,
                 counts.dim = counts.dim,
                 design.mat = x,
                 lib.sizes = lib.sizes,
                 seed = seed,
                 v.pvals = v.pvals,
                 p.pvals.1s = p.pvals.1s,
                 mu.hat.mat0 =  mu.hat.mat0,
                 phi.hat.mat0 = phi.hat.mat0,
                 ord.res.sim.mat = ord.res.sim.mat,
                 dist.mat = dist.mat,
                 pear.mat = pear.mat,
                 sim = sim
                 )
  
  # save the object as a "gofm" class
  class(gof.obj) = "gofm"
  return(gof.obj)
}
