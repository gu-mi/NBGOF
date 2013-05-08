\name{model_edgeR_trended}
\alias{model_edgeR_trended}
\title{Modeling NB trended dispersion model with the adjusted profile likelihood
estimator (APLE) on original and simulated datasets}
\usage{
  model_edgeR_trended(counts, x, lib.sizes=colSums(counts))
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{min.n}{minimim number of genes in a bin. Default is
  100. See \code{\link{dispBinTrend}} for details.}
}
\value{
  A list of quantities to be used in the main
  \code{\link{nb_gof_m}} function.
}
\description{
  This function is designed to fit an NB regression model
  with trended dispersions using the adjusted profile
  likelihood estimator. In edgeR, this function assumes the
  dispersion \eqn{\phi_i} satisfies
  \eqn{\phi_i=s(\bar{\mu}_{i\cdot})}, where \eqn{s(\cdot)}
  is a smooth function of each gene's average read counts
  across samples. A variety of non-parametric approaches
  can be used by fitting loess or spline curves on binned
  genes, or using locally weighted APL. See details below.
  The output of this function will be passed to the main
  GOF function \code{\link{nb_gof_m}}.
}
\details{
  Under the NB model ... See the
  \code{\link{eestimateGLMTrendedDisp}} and
  \code{\link{glmFit}} functions in the \code{\link{edgeR}}
  package for more information.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}
\seealso{
  \code{\link{model_edgeR_common}} and
  \code{\link{model_edgeR_genewise}}
}

