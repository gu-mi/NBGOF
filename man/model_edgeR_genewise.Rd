\name{model_edgeR_genewise}
\alias{model_edgeR_genewise}
\title{Modeling NB genewise dispersion model with the adjusted profile likelihood
estimator (APLE) on original and simulated datasets}
\usage{
  model_edgeR_genewise(counts, x,
    lib.sizes=colSums(counts))
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
  likelihood estimator. This function differs from the
  \code{\link{model_edgeR_tagwise}} function: the tagwise
  model shrinks the dispersion towards a common dispersion
  or a global dispersion trend, while the genewise model
  (implemented by this function) applies no shrinkage. See
  details below. The output of this function will be passed
  to the main GOF function \code{\link{nb_gof_m}}.
}
\details{
  In this genewise dispersion model, the dispersion
  parameter $\phi_i$ is estimated by maximizing a penalized
  log-likelihood $APL_g(\phi_g)$ without the weighted
  shared log-likelihood used in the tagwise dispersion
  model. The weight denoted by $G_0$ controls the level of
  shrinkage of purely genewise dispersions towards a common
  dispersion. In the paper of McCarthy $et al.$ (2012),
  $G_0=20/df$, where the numerator 20 is the prior degrees
  of freedom, and the denominator "df" is the residual
  degrees of freedom (number of libraries minus the number
  of distinct treatment groups). We set the value of
  \code{prior.df} equals 0, so that no shrinkage is
  applied. See the \code{\link{estimateGLMCommonDisp}},
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
  \code{\link{model_edgeR_tagwise}} for the shrinkage
  version
}

