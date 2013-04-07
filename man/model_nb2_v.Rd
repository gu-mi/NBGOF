\name{model_nb2_v}
\alias{model_nb2_v}
\title{Modeling NB2 regression model with MLE on original and simulated
datasets}
\usage{
  model_nb2_v(y, x, lib.sizes=NULL)
}
\arguments{
  \item{y}{an n-by-1 vector of non-negative integers. For a
  typical RNA-Seq experiment, this may represent the read
  counts for a single gene}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is 1 for all samples}
}
\value{
  A list of quantities to be used in the main
  \code{\link{nb_gof_v}} function.
}
\description{
  This function is designed to fit an NB2 regression model.
  The output of this function will be passed to the main
  GOF function.
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

