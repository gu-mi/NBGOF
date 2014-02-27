\name{model.edgeR.genewise2}
\alias{model.edgeR.genewise2}
\title{Modeling NB Genewise Dispersion with the Adjusted Profile Likelihood
Estimator (APLE) on Original and Simulated Datasets}
\usage{
  model.edgeR.genewise2(counts, x,
    lib.sizes=colSums(counts), min.n=min.n, method=method)
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of an RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{min.n}{minimim number of genes in a bin. Default is
  100. See \code{\link{dispBinTrend}} for details
  (lower-level function of
  \code{\link{estimateGLMTrendedDisp}}).}

  \item{method}{method for estimating the trended
  dispersion, including "auto", "bin.spline", "bin.loess",
  "power" and "spline". If NULL, then the "auto" method.
  See \code{\link{estimateGLMTrendedDisp}} for more
  details.}
}
\value{
  A list of quantities to be used in the main
  \code{\link{nb.gof.m}} function.
}
\description{
  This function fits an NB regression model with genewise
  dispersions using the adjusted profile likelihood
  estimator. This function differs from the
  \code{\link{model.edgeR.tagcom}} and
  \code{\link{model.edgeR.tagtrd}} functions: the tagwise
  model shrinks the dispersion towards a common dispersion
  or a global dispersion trend, while the genewise model
  (implemented by this function) applies \strong{no
  shrinkage}. This model is \emph{not} recommended, and is
  used for diagnostics only. See details below. The output
  of this function will be passed to the main GOF function
  \code{\link{nb.gof.m}}.
}
\details{
  In this genewise dispersion model, the dispersion
  parameter \eqn{\phi_i} is estimated by maximizing a
  penalized log-likelihood \eqn{APL_g(\phi_g)} plus a
  weighted shared log-likelihood. The weight denoted by
  \eqn{G_0} controls the level of shrinkage applied to
  purely genewise dispersions towards a common dispersion.
  In the paper of McCarthy et al. (2012), \eqn{G_0=20/df}
  performs well over real RNA-Seq datasets. Here the
  numerator 20 is the prior degrees of freedom, which can
  be specified by the argument \code{prior.df} in the
  lower-level function
  \code{\link{dispCoxReidInterpolateTagwise}} in
  \code{edgeR}. The denominator "df" is the residual
  degrees of freedom (number of libraries minus the number
  of distinct treatment groups). We set the value of
  \code{prior.df} equals 0, so that no shrinkage is
  applied. See the \code{\link{estimateGLMCommonDisp}},
  \code{\link{estimateGLMTagwiseDisp}},
  \code{\link{dispCoxReidInterpolateTagwise}} and
  \code{\link{glmFit}} functions in the \code{edgeR}
  package for more details.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  See \url{https://github.com/gu-mi/NBGOF/wiki/} for more
  details.
}
\seealso{
  \code{\link{model.edgeR.tagcom}} and
  \code{\link{model.edgeR.tagtrd}} for the shrinkage
  versions.
}

