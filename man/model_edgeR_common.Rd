\name{model_edgeR_common}
\alias{model_edgeR_common}
\title{Modeling NB2 Common Dispersion with the Adjusted Profile Likelihood
Estimator (APLE) on Original and Simulated Datasets}
\usage{
  model_edgeR_common(counts, x, lib.sizes=colSums(counts),
    design=design)
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of an RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

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
  This function fits an NB regression model with a common
  dispersion using the adjusted profile likelihood
  estimator. See details below. The output of this function
  will be passed to the main GOF function
  \code{\link{nb_gof_m}}.
}
\details{
  This function calls the
  \code{\link{estimateGLMCommonDisp}} function in
  \code{edgeR}, with the default estimation method using
  Cox-Reid adjusted profile likelihood
  (\code{\link{dispCoxReid}}). See the
  \code{\link{estimateGLMCommonDisp}} and
  \code{\link{glmFit}} functions in the \code{\link{edgeR}}
  package for more information.
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
  \code{\link{model_edgeR_tagcom}} and
  \code{\link{model_edgeR_genewise}}
}

