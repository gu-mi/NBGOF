\name{model_edgeR_tagwise}
\alias{model_edgeR_tagwise}
\title{Modeling NB tagwise dispersion model with the adjusted profile likelihood
estimator (APLE) on original and simulated datasets}
\usage{
  model_edgeR_tagwise(counts, x, lib.sizes=colSums(counts))
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
  This function is designed to fit an NB regression model
  with genewise dispersions using the adjusted profile
  likelihood estimator. In edgeR, this tagwise dispersion
  model actually lie between two extreme scenarios:
  entirely individual \eqn{\phi_i} for each gene and
  entirely shared values. This function shrinkages the
  dispersion towards the global dispersion trend. See
  details below. The output of this function will be passed
  to the main GOF function.
}
\details{
  Under the NB model ... See the
  \code{\link{estimateGLMCommonDisp}},
  \code{\link{estimateGLMTagwiseDisp}} and
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
  \code{\link{model_edgeR_genewise}} for the non-shrinkage
  version
}

