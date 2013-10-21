################################################################################
## Copyright (C) 2013 Gu Mi <mig@stat.oregonstate.edu>
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA
################################################################################

# Package Documentation
# 
# Author: Gu Mi
################################################################################

#' @title Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of 
#' RNA Sequencing Data
#' 
#' @description An R package for implementing goodness-of-fit tests for negative binomial (NB)
#' distributions and NB dispersion models, with applications in RNA-Seq data analysis.
#' 
#' @details 
#' This package can be used to test the goodness-of-fit of NB2, NBP or Poisson regression
#' models. It can also be used to test the goodness-of-fit for a variety of negative binomial
#' dispersion models in popular R/Bioconductor packages, including
#' \itemize{
#' \item NBP dispersion model in the \code{NBPSeq} package (\code{NBP})
#' \item NBQ dispersion model in the \code{NBPSeq} package (\code{NBQ})
#' \item NB common dispersion model in the \code{\link{edgeR}} package (\code{edgeR-common})
#' \item NB genewise dispersion model in the \code{\link{edgeR}} package (\code{edgeR-genewise})
#' \item NB trended (non-parametric) dispersion model in the \code{\link{edgeR}} package
#' (\code{edgeR-trended})
#' \item NB tagwise-common dispersion model in the \code{\link{edgeR}} package 
#' (\code{edgeR-tagcom})
#' \item NB tagwise-trended dispersion model in the \code{\link{edgeR}} package 
#' (\code{edgeR-tagtrd})
#' }
#' 
#' @name NBGOF-package
#' @aliases NBGOF-package NBGOF
#' @docType package
#' @import NBPSeq edgeR numDeriv parallel foreach
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' Maintainer: Gu Mi <http://people.oregonstate.edu/~mig>
#' 
#' @seealso See \code{\link{nb.gof.v}} and \code{\link{nb.gof.m}} for examples on simulated
#' datasets. See \code{\link{earthquake}} and \code{\link{arab}} for real data examples.
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' @keywords package 
#' 
NULL

.onLoad <- function(libname, pkgname){
  
  # suppress loading package messages
  suppressMessages(library(doMC))
  suppressMessages(library(foreach))
  suppressMessages(library(iterators))
  suppressMessages(library(parallel))
  #
  # startup message of our own package
#   packageStartupMessage("Loading NBGOF ...", appendLF = FALSE)
#   Sys.sleep(1)
#   packageStartupMessage(" done")
  message("\n For updates of the NBGOF package, please visit https://github.com/gu-mi/NBGOF \n")
}  

