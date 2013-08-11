
#' @title Modeling NBP Genewise Dispersion with the Maximum Likelihood Estimator 
#' (MLE) on Original and Simulated Datasets
#' 
#' @description This function fits an NBP dispersion model where the dispersion parameter
#' is modeled as a linear function of the relative means. See details below. 
#' The output of this function will be passed to the main GOF function \code{\link{nb_gof_m}}.
#' 
#' @details Under the NB model, the mean-variance relationship of a single read count 
#' satisfies \eqn{\sigma_{ij}^2 = \mu_{ij} + \phi_{ij} \mu_{ij}^2}. For applying the NBP 
#' model to RNA-Seq data, we consider the "log-linear-rel-mean" method assuming a 
#' parametric dispersion model \eqn{\phi_{ij} = \alpha_0 + \alpha_1 \log(\pi_{ij})},
#' where \eqn{\pi_{ij} = \mu_{ij}/(N_j R_j)} is the relative mean frequency after 
#' normalization. The parameters \eqn{(\alpha_0, \alpha_1)} in this dispersion model 
#' are estimated by maximizing the adjusted profile likelihood. See the 
#' \code{\link{estimate.dispersion}} function in the \code{NBPSeq} package
#' for more information.
#' 
#' @usage
#' model_nbp_m(counts, x, lib.sizes=colSums(counts))
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param lib.sizes library sizes of an RNA-Seq experiment. Default is the column
#' sums of the \code{counts} matrix.
#' 
#' @return A list of quantities to be used in the main \code{\link{nb_gof_m}} function.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references 
#' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
#' Model for Assessing Differential Gene Expression from RNA-Seq", \emph{Statistical 
#' Applications in Genetics and Molecular Biology}, 10 (1).
#'  
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
####################################################################
# model_nbp_m= function(counts, x, lib.sizes=colSums(counts)){
model_nbp_m_1 <- function(counts, x, lib.sizes=colSums(counts), obs.fit=obs.fit){
######################################################################
  nc = dim(counts)[2]
  
  # preconditions
  stopifnot(is.matrix(x), nc == dim(x)[1])
  
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  
  # data preparations
  nb.data = prepare.nb.data(counts, lib.sizes=lib.sizes)
##################################################################################### 
#  fit = estimate.dispersion(nb.data, x, print.level=0)
  # phi is obtained by estimator of alpha0, alpha1 from original data
  # alpha0 and alpha1 are extracted from "obs.fit"
  fit <-  estimate.disp.mapl.nbp.1(y=nb.data$counts, lib.sizes=lib.sizes, x=x, print.level=0,
                                   obs.fit = obs.fit)
##################################################################################### 
  phi = fit$models[[1]]$phi        # NBP "phi" --> mu+phi*mu^2
  mu = fit$models[[1]]$mu
  v = mu + phi * mu^2              # variance matrix
  res.m = (counts - mu) / sqrt(v)  # res. matrix
  
  # sort res.m with care!
  res.om = t(apply(res.m, 1, sort))
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_nbp_m_obj = list(mu.hat.mat = mu,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi
  )
  return(model_nbp_m_obj)

}
