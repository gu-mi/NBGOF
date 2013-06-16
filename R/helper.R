
################################################################################
# helper functions in NBGOF, including discarded functions for future use
# NOTE: all functions should be made invisible to end-users!
################################################################################

#' @title Group-wise sort of residual matrix in multivariate case
#' @description Group-wise sort of residual matrix in multivariate case. THIS MAY BE USEFUL
#' WHEN WE CONSIDER MULTIPLE GROUPS AND COMPLEX DESIGNS: HOW WE SORT THE RESIDUALS? CURRENTLY
#' NOT USED.
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
sort.vec = function(x, grp.ids)  ave(x, grp.ids, FUN = sort)

printf = function(...) cat(sprintf(...));
  
#' @title Based on p-values obtained, perform chi-square-like tests for overall GOF
#' 
#' @description This function is designed to provide different methods of combining p-values
#' 
#' @param x an object of class \code{gofm}
#' @param ... for future use
#' 
#' @return P-values from different methods
#' 
#' @details When the response is a count matrix, we can use this function to test
#' the goodness-of-fit of a specified negative binomial dispersion model. 
#' 
#' @usage chisq_gof(x, ...)
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
chisq_gof = function(x, ...){
  
  ## extract quantities from "gofm" object
  sim = x$sim
  v.pvals = x$v.pvals  # p-values based on vertical distances
  p.pvals.1s = x$p.pvals.1s  # p-values based on Pearson statistics (1-sided)
  p.pvals.2s = x$p.pvals.2s  # p-values based on Pearson statistics (2-sided)
  ## to avoid p-value = 0 for a two-sided test (log(p) = -Inf, so Fisher fails)
#   v.pvals[v.pvals == 0] = 1/sim
#   p.pvals.1s[p.pvals.1s == 0] = 1/sim
#   p.pvals.2s[p.pvals.2s == 0] = 2/sim
  m = length(v.pvals)
  
  ## Fisher's method of combining p-values
  # vertical distance measure
  Fisher.method.X2.vert = -2*sum(log(v.pvals))
  Fisher.method.p.vert = pchisq(Fisher.method.X2.vert, 2*m, lower.tail=FALSE) 
  # Pearson statistics (1-sided)
  Fisher.method.X2.pear.1s = -2*sum(log(p.pvals.1s))
  Fisher.method.p.pear.1s = pchisq(Fisher.method.X2.pear.1s, 2*m, lower.tail=FALSE) 
  # Pearson statistics (2-sided)
  Fisher.method.X2.pear.2s = -2*sum(log(p.pvals.2s))
  Fisher.method.p.pear.2s = pchisq(Fisher.method.X2.pear.2s, 2*m, lower.tail=FALSE) 
  
  # save as a list
  results = list(fisher.vert = Fisher.method.p.vert,
                 fisher.pear.1s = Fisher.method.p.pear.1s,
                 fisher.pear.2s = Fisher.method.p.pear.2s
                 )
  return(results)
  
}

# Stouffer.test <- function(p, w) { # p is a vector of p-values
#   if (missing(w)) {
#     w <- rep(1, length(p))/length(p)
#   } else {
#     if (length(w) != length(p))
#       stop("Length of p and w must equal!")
#   }
#   Zi <- qnorm(1-p) 
#   Z  <- sum(w*Zi)/sqrt(sum(w^2))
#   p.val <- 1-pnorm(Z)
#   return(c(Z = Z, p.value = p.val))
# }


phi.line = function (mu, v, alpha = 2, all = TRUE, ...) 
{
  id = order(mu)
  mu = mu[id]
  v = v[id]
  phi = (v - mu)/mu^alpha
  #
  if (all){
    if (length(mu) > 1000) {
      id = c(seq(1, length(mu) - 1, length = 1000), length(mu))
    }
    lines(mu, phi, ...)
    invisible()
  }
  else if (!all){
    id = (phi > 0) & (mu > 0)
    mu = mu[id]
    phi = phi[id]
    if (length(mu) > 1000) {
      id = c(seq(1, length(mu) - 1, length = 1000), length(mu))
    }
    lines(mu[id], phi[id], ...)
    invisible()
  }
  
}

## scatter plot of data points ##
phi.sp = function(counts, all = TRUE, alpha = 2, ...){
  v = apply(counts, 1, var)
  mu = apply(counts, 1, mean)
  phi = (v - mu)/mu^alpha
  if (all){
    plot(mu, phi, pch=".", cex=3, ...)
  }
  else if (!all){
    id = (phi > 0) & (mu > 0)
    mu = mu[id]
    phi = phi[id]
    plot(mu, phi, pch=".", cex=3, ...)
  }
  
}

plot.mphi = function(counts, fitted.lines,
                     cols = 1:length(fitted.lines), lwds = 2, ltys = 1, alpha = 2, ...) {
  
  ## fitted.lines:  a list of fitted lines
  
  n.lines = length(fitted.lines)
  
  cols = rep(cols, n.lines)[1:n.lines]
  lwds = rep(lwds, n.lines)[1:n.lines]
  ltys = rep(ltys, n.lines)[1:n.lines]
  
  phi.sp(counts, alpha=alpha, ...)
  
  if (n.lines>0) {
    for (i in 1:n.lines) {
      phi.line(fitted.lines[[i]]$mu, fitted.lines[[i]]$v,
               col=cols[i], lty=ltys[i], lwd=lwds[i], alpha=alpha)
    }
  }
}

# combine.p --- A short and ugly function to combine probabilities
# from independent tests of significance using either the Fisher,
# Winer or Stouffer's methods
# see:
# Fisher (1958, section 21.1 --- pp 99--101)
# Sokal & Rohlf (1995, section 18.1, pp 794--797)
# Wolf (1986, pp 18--23)
#
# To do:
# 0) fix winer & stouffer to handle signed cases, this is broken!
# 1) add bulletproofing for y stuff, as it stands now, winer
#    will accept all sorts of crappy input
# 2) Beautify using switch on method, or not...
#
# Pete Hurd June 20 2002

# combine.p <- function(x, y = NULL, method="fisher") {
#   x <- as.vector(na.omit(x))
#   if (any(x <= 0) || any(is.na(x)) || any(x>1))
#     stop("all entries of x must be in the range (0,1]")
#   
#   if( method=="fisher"){
#     STATISTIC <- 2*(sum(-1*log(x)))
#     PARAMETER <- 2*length(x)
#     PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
#     METHOD <- "Fisher's combined probability"
#     names(STATISTIC) <- "Chi Square"
#     names(PARAMETER) <- "X-squared df"
#     names(PVAL) <- "p.value"
#   }
#   else if (method=="winer"){
#     t <- qt(x,y,lower.tail=FALSE)
#     STATISTIC <- sum(t) / (sqrt(sum(dv/(dv-2))))
#     PARAMETER <- NA
#     PVAL <-pnorm(STATISTIC,lower.tail=FALSE)
#     METHOD <- "Winer's combined probability"
#     names(STATISTIC) <- "Z score"
#     names(PARAMETER) <- "No df"
#     names(PVAL) <- "p.value"    
#   }
#   else if (method=="stouffer"){
#     STATISTIC <- sum(qnorm(x,lower.tail=FALSE)) / sqrt(length(x))
#     PARAMETER <- NA
#     PVAL <-pnorm(STATISTIC,lower.tail=FALSE)
#     METHOD <- "Stouffer's combined probability"
#     names(STATISTIC) <- "Z score"
#     names(PARAMETER) <- "No df"
#     names(PVAL) <- "p.value"
#   }
#   
#   structure(list(statistic=STATISTIC, parameter=PARAMETER,
#                  p.value=PVAL,method=METHOD), class="htest")
# }


# http://www.stat.wisc.edu/~st571-1/gtest.R

# g.test <- function(x, y = NULL, correct="williams",
#                    p = rep(1/length(x), length(x)), simulate.p.value = FALSE, B = 2000)
#   #can also use correct="none" or correct="yates"
# {
#   DNAME <- deparse(substitute(x))
#   if (is.data.frame(x)) x <- as.matrix(x)
#   if (is.matrix(x)) {
#     if (min(dim(x)) == 1) 
#       x <- as.vector(x)
#   }
#   if (!is.matrix(x) && !is.null(y)) {
#     if (length(x) != length(y)) 
#       stop("x and y must have the same length")
#     DNAME <- paste(DNAME, "and", deparse(substitute(y)))
#     OK <- complete.cases(x, y)
#     x <- as.factor(x[OK])
#     y <- as.factor(y[OK])
#     if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
#       stop("x and y must have at least 2 levels")
#     x <- table(x, y)
#   }
#   if (any(x < 0) || any(is.na(x))) 
#     stop("all entries of x must be nonnegative and finite")
#   if ((n <- sum(x)) == 0) 
#     stop("at least one entry of x must be positive")
#   #If x is matrix, do test of independence
#   if (is.matrix(x)) {
#     #Test of Independence
#     nrows<-nrow(x)
#     ncols<-ncol(x)
#     if (correct=="yates"){ # Do Yates' correction?
#       if(dim(x)[1]!=2 || dim(x)[2]!=2) # check for 2x2 matrix
#         stop("Yates' correction requires a 2 x 2 matrix")
#       if((x[1,1]*x[2,2])-(x[1,2]*x[2,1]) > 0)
#       {
#         x[1,1] <- x[1,1] - 0.5
#         x[2,2] <- x[2,2] - 0.5
#         x[1,2] <- x[1,2] + 0.5
#         x[2,1] <- x[2,1] + 0.5
#       }
#       else
#       {
#         x[1,1] <- x[1,1] + 0.5
#         x[2,2] <- x[2,2] + 0.5
#         x[1,2] <- x[1,2] - 0.5
#         x[2,1] <- x[2,1] - 0.5
#       }
#     }
#     
#     sr <- apply(x,1,sum)
#     sc <- apply(x,2,sum)
#     E <- outer(sr,sc, "*")/n
#     # are we doing a monte-carlo?
#     # no monte carlo GOF?
#     if (simulate.p.value){
#       METHOD <- paste("Log likelihood ratio (G-test) test of independence\n\t with simulated p-value based on", B, "replicates")
#       tmp <- .C("gtestsim", as.integer(nrows), as.integer(ncols),
#                 as.integer(sr), as.integer(sc), as.integer(n), as.integer(B),
#                 as.double(E), integer(nrows * ncols), double(n+1),
#                 integer(ncols), results=double(B), PACKAGE= "ctest")
#       g <- 0
#       for (i in 1:nrows){
#         for (j in 1:ncols){
#           if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
#         }
#       }
#       STATISTIC <- G <- 2 * g
#       PARAMETER <- NA
#       PVAL <- sum(tmp$results >= STATISTIC)/B
#     }
#     else {
#       # no monte-carlo
#       # calculate G
#       g <- 0
#       for (i in 1:nrows){
#         for (j in 1:ncols){
#           if (x[i,j] != 0) g <- g + x[i,j] * log(x[i,j]/E[i,j])
#         }
#       }
#       q <- 1
#       if (correct=="williams"){ # Do Williams' correction
#         row.tot <- col.tot <- 0    
#         for (i in 1:nrows){ row.tot <- row.tot + 1/(sum(x[i,])) }
#         for (j in 1:ncols){ col.tot <- col.tot + 1/(sum(x[,j])) }
#         q <- 1+ ((n*row.tot-1)*(n*col.tot-1))/(6*n*(ncols-1)*(nrows-1))
#       }
#       STATISTIC <- G <- 2 * g / q
#       PARAMETER <- (nrow(x)-1)*(ncol(x)-1)
#       PVAL <- 1-pchisq(STATISTIC,df=PARAMETER)
#       if(correct=="none")
#         METHOD <- "Log likelihood ratio (G-test) test of independence without correction"
#       if(correct=="williams")
#         METHOD <- "Log likelihood ratio (G-test) test of independence with Williams' correction"
#       if(correct=="yates")
#         METHOD <- "Log likelihood ratio (G-test) test of independence with Yates' correction"
#     }
#   }
#   else {
#     # x is not a matrix, so we do Goodness of Fit
#     METHOD <- "Log likelihood ratio (G-test) goodness of fit test"
#     if (length(x) == 1) 
#       stop("x must at least have 2 elements")
#     if (length(x) != length(p)) 
#       stop("x and p must have the same number of elements")
#     E <- n * p
#     
#     if (correct=="yates"){ # Do Yates' correction
#       if(length(x)!=2)
#         stop("Yates' correction requires 2 data values")
#       if ( (x[1]-E[1]) > 0.25) {
#         x[1] <- x[1]-0.5
#         x[2] <- x[2]+0.5
#       }
#       else if ( (E[1]-x[1]) > 0.25){
#         x[1] <- x[1]+0.5
#         x[2] <- x[2]-0.5
#       }
#     }
#     names(E) <- names(x)
#     g <- 0
#     for (i in 1:length(x)){
#       if (x[i] != 0) g <- g + x[i] * log(x[i]/E[i])
#     }
#     q <- 1
#     if (correct=="williams"){ # Do Williams' correction
#       q <- 1+(length(x)+1)/(6*n)
#     }
#     STATISTIC <- G <- 2*g/q
#     PARAMETER <- length(x) - 1
#     PVAL <- pchisq(STATISTIC, PARAMETER, lower = FALSE)
#   }
#   names(STATISTIC) <- "Log likelihood ratio statistic (G)"
#   names(PARAMETER) <- "X-squared df"
#   names(PVAL) <- "p.value"
#   structure(list(statistic=STATISTIC,parameter=PARAMETER,p.value=PVAL,
#                  method=METHOD,data.name=DNAME, observed=x, expected=E),
#             class="htest")
# }


#' @title User-defined link function for gamma regression
#' 
#' @description User-defined link function for gamma regression. This is made available to be passed
#' to the \code{family=} argument in the \code{glm} function
#' 
#' @usage
#' family = Gamma(link = vlog())
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @keywords internal
#' 
#' @examples
#' 
#' library(NBGOF)
#' set.seed(123)
#' n = 1000                     
#' x = runif(n)
#' sh = 2   
#' vv = vlog()                     
#' y = rgamma(n, scale=vv$linkinv(2+3*x)/sh, shape=sh)
#' fit1 = glm(y~x, family=Gamma(link=vv))                       
#' fit1
#' sum.fit1 = summary(fit1)
#' 1/sum.fit1$dispersion   # compare with sh: close to each other as n increases
#' 
vlog <- function() {
  ## link
  linkfun <- function(y) log(exp(y)-1)
  ## inverse link
  linkinv <- function(eta)  log(exp(eta)+1)
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 1/(exp(-eta) + 1) }
  valideta <- function(eta) TRUE
  link <- "log(exp(y)-1)"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}


