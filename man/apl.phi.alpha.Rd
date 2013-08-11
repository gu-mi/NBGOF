\name{apl.phi.alpha}
\alias{apl.phi.alpha}
\title{The log adjusted profile likelihood of the dispersion parameters from all genes}
\usage{
  apl.phi.alpha(phi0, alpha1, y, lib.sizes, x, beta0,
    subset = 1:m, per.row = FALSE, print.level = 0)
}
\arguments{
  \item{phi0}{a number, dispersion model parameter}

  \item{alpha1}{a number, dispsersion model parameter}

  \item{y}{an m*n matrix of counts}

  \item{lib.sizes}{an n vector of effective library sizes}

  \item{x}{a n by p design matrix}

  \item{beta0}{non NA components are hypothesized values of
  beta, NA components are free components}

  \item{subset}{specify the subset of rows to use to
  compute the likelihood}

  \item{per.row}{logical, if \code{TRUE} return the
  likelihood from each row, otherwise return the total
  likelihood}

  \item{print.level}{a number specifying print level}
}
\value{
  the log adjusted profile likelihood of phi
}
\description{
  The log adjusted profile likelihood of the dispersion
  parameters from all genes
}

