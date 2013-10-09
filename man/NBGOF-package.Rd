\docType{package}
\name{NBGOF-package}
\alias{NBGOF}
\alias{NBGOF-package}
\title{Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of
RNA Sequencing Data}
\description{
  An R package for implementing goodness-of-fit tests for
  negative binomial (NB) distributions and NB dispersion
  models, with applications in RNA-Seq data analysis.
}
\details{
  This package can be used to test the goodness-of-fit of
  NB2, NBP or Poisson regression models. It can also be
  used to test the goodness-of-fit for a variety of
  negative binomial dispersion models in popular
  R/Bioconductor packages, including \itemize{ \item NBP
  dispersion model in the \code{NBPSeq} package
  (\code{NBP}) \item NBQ dispersion model in the
  \code{NBPSeq} package (\code{NBQ}) \item NB common
  dispersion model in the \code{\link{edgeR}} package
  (\code{edgeR-common}) \item NB genewise dispersion model
  in the \code{\link{edgeR}} package
  (\code{edgeR-genewise}) \item NB trended (non-parametric)
  dispersion model in the \code{\link{edgeR}} package
  (\code{edgeR-trended}) \item NB tagwise-common dispersion
  model in the \code{\link{edgeR}} package
  (\code{edgeR-tagcom}) \item NB tagwise-trended dispersion
  model in the \code{\link{edgeR}} package
  (\code{edgeR-tagtrd}) }
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer

  Maintainer: Gu Mi <http://people.oregonstate.edu/~mig>
}
\references{
  See \url{https://github.com/gu-mi/NBGOF/wiki/} for more
  details.
}
\seealso{
  See \code{\link{nb.gof.v}} and \code{\link{nb.gof.m}} for
  examples on simulated datasets. See
  \code{\link{earthquake}} and \code{\link{arab}} for real
  data examples.
}
\keyword{package}

