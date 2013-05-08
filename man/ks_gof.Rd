\name{ks_gof}
\alias{ks_gof}
\title{Main Function of Implementing Simulation-based Goodness-of-Fit Tests on
the regression of a univariate count response on one or more explanatory variables}
\usage{
  ks_gof(y0, x, model = "NB2", sim=1)
}
\arguments{
  \item{y0}{an n-by-1 vector of non-negative integers. For
  a typical RNA-Seq experiment, this may represent the read
  counts for a single gene}

  \item{x}{an n-by-p design matrix.}

  \item{sim}{number of simulations performed.}

  \item{model}{a string of characters specifying the
  negative binomial model used to fit the data. Currently
  the following dispersion models are available to be
  checked for goodness-of-fit}
}
\value{
  An object
}
\description{
  This function is designed to test goodness-of-fit of an
  NB2 or NBP negative binomial regression model with a
  known design matrix. Estimation methods for NBP model
  fitting include MLE and APLE
}
\details{
  When the response is a vector, we can use this function
  to test the goodness-of-fit of a specified negative
  binomial regression model.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}

