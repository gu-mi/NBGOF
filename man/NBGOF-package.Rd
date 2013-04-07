\docType{package}
\name{NBGOF-package}
\alias{NBGOF}
\alias{NBGOF-package}
\title{Simulation-based Goodness-of-Fit Tests for Evaluating the Adequacy of
Negative Binomial Regression/Dispersion Models}
\description{
  Motivated by the needs for evaluating negative binomial
  (NB) model adequacy based on different dispersion models,
  this package provides a simulation-based approach to
  check the goodness-of-fit for a given NB regression
  model.
}
\details{
  It can be used to test the goodness-of-fit for a variety
  of negative binomial dispersion models in popular
  R/Bioconductor packages, including \itemize{ \item NB
  genewise dispersion model in the \code{\link{NBPSeq}}
  package \item NBP dispersion model in the
  \code{\link{NBPSeq}} package \item NB common dispersion
  model in the \code{\link{edgeR}} package \item NB tagwise
  dispersion model in the \code{\link{edgeR}} package \item
  NB trended dispersion model in the \code{\link{edgeR}}
  package }
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer

  Maintainer: Gu Mi <http://people.oregonstate.edu/~mig>
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}
\seealso{
  See \code{\link{nb_gof_v}} and \code{\link{nb_gof_m}} for
  examples on simulated datasets
}
\keyword{package}

