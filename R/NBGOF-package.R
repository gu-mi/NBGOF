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

#' Simulation-based Goodness-of-Fit Tests for Evaluating the Adequacy of
#' Negative Binomial Regression/Dispersion Models
#' 
#' Motivated by the needs for evaluating negative binomial (NB) model adequacy based on
#' different dispersion models, this package provides a simulation-based approach to 
#' check the goodness-of-fit for a given NB regression model.
#' 
#' It can be used to test the goodness-of-fit for a variety of negative binomial
#' dispersion models in popular R/Bioconductor packages, including
#' \itemize{
#' \item NB genewise dispersion model in the \code{\link{NBPSeq}} package
#' \item NBP dispersion model in the \code{\link{NBPSeq}} package
#' \item NB common dispersion model in the \code{\link{edgeR}} package
#' \item NB tagwise dispersion model in the \code{\link{edgeR}} package
#' \item NB trended dispersion model in the \code{\link{edgeR}} package
#' }
#' 
#' @name NBGOF-package
#' @aliases NBGOF-package NBGOF
#' @docType package
#' @import NBPSeq edgeR
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' Maintainer: Gu Mi <http://people.oregonstate.edu/~mig>
#' 
#' @seealso See \code{\link{nb_gof_v}} and \code{\link{nb_gof_m}} for examples on simulated
#' datasets
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' @keywords package 
#' 
NULL

.onLoad <- function(libname, pkgname){
  message <- "\n The NBGOF package is experimental and under development: 
  Function syntax may have minor changes in future versions.
  See https://github.com/gu-mi/NBGOF for more information. \n"
  packageStartupMessage(message)
}  

