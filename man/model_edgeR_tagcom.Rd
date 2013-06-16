\name{model_edgeR_tagcom}
\alias{model_edgeR_tagcom}
\title{Modeling NB Tagwise-Common Dispersion with the Adjusted Profile Likelihood
Estimator (APLE) on Original and Simulated Datasets}
\usage{
  model_edgeR_tagcom(counts, x, lib.sizes=colSums(counts),
    prior.df = prior.df, design=design)
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of an RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{prior.df}{prior degrees of freedom to control the
  level of shrinkage (default is 10). This argument is
  controlled by a higher level function
  \code{\link{nb_gof_m}}.}

  \item{design}{specifications of testing (1) a
  single-group model (\code{single}); (2) a multiple-group
  model (\code{multiple}); (3) complex design with
  interactions (\code{complex}). \strong{The codes are only
  tested for single- and multiple-group cases, and the
  modeling of dispersions in the multiple-group case is
  still under consideration! For the complex design with
  interactions, though we can pass the design matrix
  directly in edgeR, the results may not make any sense!}}
}
\value{
  A list of quantities to be used in the main
  \code{\link{nb_gof_m}} function.
}
\description{
  This function fits an NB regression model with tagwise
  dispersions using the adjusted profile likelihood
  estimator (APLE). In \code{edgeR}, this tagwise
  dispersion model lies between two extreme scenarios:
  entirely individual \eqn{\phi_i} for each gene and
  entirely shared values (i.e. common dispersion). The
  function \code{estimateGLMTagwiseDisp} in \code{edgeR}
  shrinks the dispersion towards a common dispersion or a
  global dispersion trend. In our implementation, the
  function \code{model_edgeR_tagcom} shrinks the tagwise
  dispersions towards a common dispersion. See details
  below. The output of this function will be passed to the
  main GOF function \code{\link{nb_gof_m}}.
}
\details{
  In this tagwise-common dispersion model, the dispersion
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
  of distinct treatment groups). We leave the default value
  of \code{prior.df} at 10. See the
  \code{\link{estimateGLMCommonDisp}},
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
  \code{\link{model_edgeR_genewise}} for the non-shrinkage,
  entirely genewise model, and
  \code{\link{model_edgeR_tagtrd}} for the tagwise model
  that shrinks towards the trended dispersion.
}

