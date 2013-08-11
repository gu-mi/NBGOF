\name{estimate.disp.mapl.nbp.1}
\alias{estimate.disp.mapl.nbp.1}
\title{Estimate the NBP model for Dispersion Parameters}
\usage{
  estimate.disp.mapl.nbp.1(y, lib.sizes, x,
    beta = rep(NA, dim(x)[2]), alpha1 = NA, mu.lower = 1,
    mu.upper = Inf,
    subset = rowSums(is.na(mu.pre) | mu.pre < mu.lower | mu.pre > mu.upper) == 0,
    mu.pre = irls.nb(y, lib.sizes, x, phi = phi.pre, beta)$mu,
    phi.pre = 0.1, tol.alpha0 = 1e-04, tol.alpha1 = 0.001,
    print.level = 1, obs.fit = obs.fit)
}
\arguments{
  \item{y}{an m*n matrix of NB counts}

  \item{lib.sizes}{effective library sizes (expected column
  totals)}

  \item{x}{an n*p model matrix}

  \item{beta}{a p-vector of coefficients in a log linear
  model for the mean counts: log(mu) = x' beta.}

  \item{alpha1}{if given}

  \item{mu.lower}{a number, rows with mu.pre < mu.lower
  will not be used for estimating the dispersion model}

  \item{mu.upper}{a number, rows with mu.pre > mu.upper
  will not be used for estimating the dispersion model}

  \item{subset}{specify the subest of rows to use to
  compute the likelihood}

  \item{mu.pre}{a matrix, preliminary estimates of the mu
  matrix}

  \item{phi.pre}{a number, a preliminary estimate of phi}

  \item{tol.alpha0}{a number specifying the convergence
  tolerance for alpha0}

  \item{tol.alpha1}{a number specifying the convergence
  tolerance for alpha1}

  \item{print.level}{a number, print level}
}
\value{
  a list
}
\description{
  Estimate the NBP model for dispersion parameters by
  maximizing the adjusted profile likelihood.
}
\details{
  Under the NBP model, the dispersion is a power function
  of the relative mean,

  phi[i,j] = phi0 (mu[i,j]/s[j])^alpha1.

  This function estiamtes the two parameters (phi0, alpha1)
  of the NBP model by maximizing the adjusted profile
  likelihood of (phi_0, alpha_1).
}
\author{
  Yanming Di
}

