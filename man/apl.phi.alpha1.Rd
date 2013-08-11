\name{apl.phi.alpha1}
\alias{apl.phi.alpha1}
\title{The log adjusted profile likelihood of the dispersion model parameters}
\usage{
  apl.phi.alpha1(phi0, alpha1, y, lib.sizes = lib.sizes, x,
    beta0, print.level = 0)
}
\arguments{
  \item{phi0}{a number, dispersion model parameter}

  \item{alpha1}{a number, dispsersion model parameter}

  \item{y}{an n-vector of NB counts}

  \item{lib.sizes}{an n vector of effective library sizes}

  \item{x}{a n by p design matrix}

  \item{beta0}{non NA components are hypothesized values of
  beta, NA components are free components}

  \item{print.level}{a number, print level}
}
\value{
  a list the log adjusted profile likelihood of (phi0,
  alpha1)
}
\description{
  The log adjusted profile likelihood of the dispersion
  model parameters from a single gene
}
\details{
  The dispersion model is phi=phi0 pi^alpha1, where pi
  denotes the mean relative frequencies
}

