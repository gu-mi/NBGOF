\name{model_nbp_m}
\alias{model_nbp_m}
\title{Modeling NBP genewise dispersion model with the maximum likelihood estimator
(MLE) on original and simulated datasets}
\usage{
  model_nbp_m(counts, x, lib.sizes=colSums(counts))
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}
}
\value{
  A list of quantities to be used in the main
  \code{\link{nb_gof_m}} function.
}
\description{
  This function is designed to fit an NBP regression model
  with genewise dispersions using the adjusted profile
  likelihood estimator. See details below. The output of
  this function will be passed to the main GOF function.
}
\details{
  Under the NB model, the mean-variance relationship of a
  single read count satisfies \eqn{\sigma_{ij}^2 = \mu_{ij}
  + \phi_{ij} \mu_{ij}^2}. For applying the NBP model to
  RNA-Seq data, we consider the "log-linear-rel-mean"
  method assuming a parametric dispersion model
  \eqn{\phi_{ij} = \alpha_0 + \alpha_1 \log(\pi_{ij})},
  where \eqn{\pi_{ij} = \mu_{ij}/(N_j R_j)} is the relative
  mean frequency after normalization. The parameters
  \eqn{(\alpha_0, \alpha_1)} in this dispersion model are
  estimated by maximizing the adjusted profile likelihood.
  See the \code{\link{estimate.dispersion}} function in the
  \code{\link{NBPSeq}} package for more information.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}

