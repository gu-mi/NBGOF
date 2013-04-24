\name{irls.nb}
\alias{irls.nb}
\title{Estiamte the regression coefficients in an NB GLM model}
\usage{
  irls.nb(y, s, x, phi, beta0, ..., print.level = 0)
}
\arguments{
  \item{y}{an m*n matrix of counts}

  \item{s}{an n vector of effective library sizes}

  \item{x}{an n*p design matrix}

  \item{phi}{an n*p matrix of regression coefficients}

  \item{tol.mu}{convergence criteria}

  \item{beta0}{a K vector, non NA components are
  hypothesized values of beta, NA components are free
  components}
}
\value{
  beta a K vector, the MLE of the regression coefficients.
}
\description{
  Estimate the regression coefficients in an NBP GLM model
  for each gene
}
\keyword{internal}

