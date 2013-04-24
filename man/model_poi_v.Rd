\name{model_poi_v}
\alias{model_poi_v}
\title{Modeling Poisson regression model with MLE on original and simulated
datasets}
\usage{
  model_poi_v(y, x, lib.sizes=NULL)
}
\arguments{
  \item{y}{an n-by-1 vector of non-negative integers. For a
  typical RNA-Seq experiment, this may represent the read
  counts for a single gene}

  \item{x}{an n-by-p design matrix. For Poisson model
  fitting, we used the \link{glm} function, so if an
  intercept is desired, there is no need to include the
  first column of 1.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is 0 for all samples}
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

