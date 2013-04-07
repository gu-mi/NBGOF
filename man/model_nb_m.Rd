\name{model_nb_m}
\alias{model_nb_m}
\title{Modeling NB genewise dispersion model with MLE on original and simulated
datasets}
\usage{
  model_nb_m(counts, x, lib.sizes=colSums(counts))
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
  with genewise dispersions using the maximum likelihood
  estimator. The output of this function will be passed to
  the main GOF function.
}
\details{
  Details here
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}

