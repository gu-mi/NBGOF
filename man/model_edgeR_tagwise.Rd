\name{model_edgeR_tagwise}
\alias{model_edgeR_tagwise}
\title{Modeling NB tagwise dispersion model with the adjusted profile likelihood
estimator (APLE) on original and simulated datasets}
\usage{
  model_edgeR_tagwise(counts, x, lib.sizes=colSums(counts),
    prior.df = prior.df)
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{prior.df}{prior degrees of freedom to control the
  level of shrinkage (default is 10)}
}
\value{
  A list of quantities to be used in the main
  \code{\link{nb_gof_m}} function.
}
\description{
  This function is designed to fit an NB regression model
  with genewise dispersions using the adjusted profile
  likelihood estimator. In \code{edgeR}, this tagwise
  dispersion model lies between two extreme scenarios:
  entirely individual \eqn{\phi_i} for each gene and
  entirely shared values (i.e. common dispersion). This
  function shrinks the dispersion towards a common
  dispersion or a global dispersion trend. See details
  below. The output of this function will be passed to the
  main GOF function \code{\link{nb_gof_m}}.
}
\details{
  In this tagwise dispersion model, the dispersion
  parameter $\phi_i$ is estimated by maximizing a penalized
  log-likelihood $APL_g(\phi_g)$ plus a weighted shared
  log-likelihood. The weight denoted by $G_0$ controls the
  level of shrinkage of purely genewise dispersions towards
  a common dispersion. In the paper of McCarthy $et al.$
  (2012), $G_0=20/df$ performs well over real RNA-Seq
  datasets. Here the numerator 20 is the prior degrees of
  freedom that can be specified (\code{prior.df}) in the
  lower-level function
  \code{\link{dispCoxReidInterpolateTagwise}} in
  \code{\link{edgeR}}. The denominator "df" is the residual
  degrees of freedom (number of libraries minus the number
  of distinct treatment groups). We leave the default value
  of \code{prior.df} at 10. See the
  \code{\link{estimateGLMCommonDisp}},
  \code{\link{estimateGLMTagwiseDisp}},
  \code{\link{dispCoxReidInterpolateTagwise}} and
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
  \code{\link{model_edgeR_genewise}} for the non-shrinkage,
  entirely genewise version
}

