\name{model.edgeR.common}
\alias{model.edgeR.common}
\title{Modeling NB2 Common Dispersion with the Adjusted Profile Likelihood
Estimator (APLE) on Original and Simulated Datasets}
\usage{
model.edgeR.common(counts, x, lib.sizes=colSums(counts), method=method)
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of an RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{method}{method for estimating the common
  dispersion, including "CoxReid", "Pearson" and
  "deviance". If NULL, then use the "CoxReid" method. See
  \code{\link{estimateGLMCommonDisp}} for more details.}
}
\value{
A list of quantities to be used in the main
\code{\link{nb.gof.m}} function.
}
\description{
This function fits an NB regression model with a common
dispersion using the adjusted profile likelihood estimator.
See details below. The output of this function will be
passed to the main GOF function \code{\link{nb.gof.m}}.
}
\details{
This function calls the \code{\link{estimateGLMCommonDisp}}
function in \code{edgeR}, with the default estimation
method using Cox-Reid adjusted profile likelihood
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
\code{\link{model.edgeR.tagcom}} and
\code{\link{model.edgeR.genewise}}
}

