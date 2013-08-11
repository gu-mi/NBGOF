## List of main functions
## 
## Likelihood (profile likelihood, conditonal likelihood, and adjusted
## profile likelihood) functions:
##
##  apl.phi.alpha.1
##  apl.phi.alpha
##
## Genewise dispersion parameter estimators:
##
##
## Estimators for dispersion parameter models
##
##  estimate.disp.mapl.nbp

##' The log adjusted profile likelihood of the dispersion model parameters from a single gene
##'
##' The dispersion model is phi=phi0 pi^alpha1, where pi denotes the mean relative frequencies
##' 
##' @title The log adjusted profile likelihood of the dispersion model parameters
##' @param phi0 a number, dispersion model parameter
##' @param alpha1 a number, dispsersion model parameter
##' @param y an n-vector of NB counts
##' @param lib.sizes an n vector of effective library sizes
##' @param x a n by p design matrix
##' @param beta0 non NA components are hypothesized values of beta, NA components are free components
##' @param print.level a number, print level
##' @return a list the log adjusted profile likelihood of (phi0, alpha1)
apl.phi.alpha.1 = function(phi0, alpha1, y, lib.sizes, x, beta0,
  print.level=0) {
  ## Estimate the regresssion coefficients in the NB regression model
  fit = irls.nbp.1(y, lib.sizes, x, phi0, alpha1, beta0);

  kappa = 1/fit$phi;

  ## Compute the information matrices
  mu.hat = fit$mu;
  v.hat = drop(mu.hat + fit$phi*mu.hat^2)
  ## i.hat=t(x)%*% diag(mu.hat^2/v.hat)%*% x;
  j.hat = t(x)%*% diag(mu.hat^2 * (y/mu.hat^2 - (y+kappa)/(mu.hat+kappa)^2) - (y-mu.hat)*mu.hat/v.hat) %*% x
  ## The likelihood of (kappa, mu) (different from the l.hat in
  ## irls.nb.1, which is the likelihood of mu)
  l.hat = log.likelihood.nb(kappa, fit$mu, y)
  ## l.hat = log.likelihood.nb.2(kappa, fit$mu, y)

  if(print.level>2) {
    print(list(mu.hat=mu.hat, v.hat=v.hat, l.hat=l.hat,
      ## i.hat=i.hat
      j.hat=j.hat))
  }

  fit$likelihood=l.hat-0.5*log(det(j.hat));

  fit
 }

 ##' @title The log adjusted profile likelihood of the dispersion parameters from all genes
##'
##' @param phi0 a number, dispersion model parameter
##' @param alpha1 a number, dispsersion model parameter
##' @param y an m*n matrix of counts
##' @param lib.sizes an n vector of effective library sizes
##' @param x a n by p design matrix
##' @param beta0 non NA components are hypothesized values of beta, NA components are free components
##' @param subset specify the subset of rows to use to compute the likelihood
##' @param per.row logical, if \code{TRUE} return the likelihood from each row, otherwise return the total likelihood 
##' @param print.level a number specifying print level
##' @return the log adjusted profile likelihood of phi

apl.phi.alpha = function(phi0, alpha1,
  y, lib.sizes, x, beta0,
  subset=1:m,
  per.row = FALSE,
  print.level=0) {

  m = dim(y)[1];
  n = dim(y)[2];

  if (length(lib.sizes)==1) {
    lib.sizes = rep(lib.sizes, n);
  }
  
  subset=(1:m)[subset];

  if (print.level>1) pb=txtProgressBar(style=3);

  l = rep(NA, m);
  ## mu = matrix(0, m, n);
  
  for (i in subset) {
    if (print.level>1) setTxtProgressBar(pb, i/m);

    obj = apl.phi.alpha.1(phi0, alpha1, y[i,], lib.sizes=lib.sizes, x, beta0,
      print.level=print.level-1);
    l[i] = obj$likelihood;
    ## mu[i,] = obj$mu;
  }
  if (print.level>1) close(pb);

  if (print.level>2)
    message(sprintf("%d genes are used for computing likelihood.", sum(subset)));

  if (per.row) l else sum(l[subset]);
}


  
##' Estimate the NBP model for dispersion parameters by maximizing the
##' adjusted profile likelihood. 
##'
##' Under the NBP model, the dispersion is a power function of the relative mean, 
##'
##'   phi[i,j] = phi0 (mu[i,j]/s[j])^alpha1.
##'
##' This function estiamtes the two parameters (phi0, alpha1) of the
##' NBP model by maximizing the adjusted profile likelihood of (phi_0,
##' alpha_1).
##'
##' @title Estimate the NBP model for Dispersion Parameters
##'
##' @param y an m*n matrix of NB counts
##' @param lib.sizes effective library sizes (expected column totals)
##' @param x an n*p model matrix
##' @param beta a p-vector of coefficients in a log linear model for
##' the mean counts: log(mu) = x' beta.
##' @param alpha1 if given 
##' @param mu.lower  a number, rows with mu.pre < mu.lower will not be used for estimating the dispersion model
##' @param mu.upper  a number, rows with mu.pre > mu.upper will not be used for estimating the dispersion model
##' @param subset specify the subest of rows to use to compute the likelihood
##' @param mu.pre a matrix, preliminary estimates of the mu matrix
##' @param phi.pre  a number, a preliminary estimate of phi
##' @param tol.alpha0 a number specifying the convergence tolerance for alpha0
##' @param tol.alpha1 a number specifying the convergence tolerance for alpha1
##' @param print.level a number, print level
##' @return a list 
##' @author Yanming Di

##############################################################
 estimate.disp.mapl.nbp.1 <- function(y,lib.sizes, x,
 # estimate.disp.mapl.nbp = function(y, lib.sizes, x,
#########################################################
  beta=rep(NA, dim(x)[2]),
  alpha1=NA,
  mu.lower=1, mu.upper=Inf,
  subset = rowSums(is.na(mu.pre) | mu.pre<mu.lower | mu.pre>mu.upper)==0,
  mu.pre = irls.nb(y, lib.sizes, x, phi=phi.pre, beta)$mu,
  phi.pre = 0.1,
  tol.alpha0 = 1e-4,
  tol.alpha1 = 1e-3,
  print.level=1,
  obs.fit=obs.fit) {  # ADD obs.fit FOR THE ORIGINIAL FIT, TO GET alpha0 and alpha1

  ## Mean relative counts
  m = dim(y)[1];
  n = dim(y)[2];
  mu.pre = matrix(mu.pre, m, n);

  mu.lower = max(mu.lower, 1);

  p = mu.pre / (matrix(1, m, 1) %*% matrix(lib.sizes, 1, n));

  if (print.level>0) {
    message("Estimating the dispersion model:");
    message("log(dispersion) = alpha0 + alpha1 log(mu/lib.size)");
  }

############################################################################
  ## Bounds for alpha0 and alpha1.
  ## We are mainly interested in alpha1 in [1,2].
  # alpha.bounds = c(1 - tol.alpha, 2 + tol.alpha)
## USE alpha0, alpha1 from estimate.dispersion in real data
  # alpha1.bounds = c(-1.5, 0.5);
  # alpha0.lower = -20;

#  N = max(lib.sizes);

 # l.alpha0 = function(alpha0, alpha1) {
 # if (print.level>2) {
 #	message(sprintf("alpha0=%f, phi=%f when mu=1000 in a library of size %.0f.\n",
 #            alpha0, exp(alpha0)*(1000/N)^alpha1, N));
 # }

#    apl.phi.alpha(exp(alpha0), alpha1, y, lib.sizes, x, beta,
#                  subset=subset,
#                  print.level=print.level-1);
#  }

   ## The adjusted profile likelihood of alpha
#  pl.alpha1 = function(alpha1) {
#    if (print.level>0) message(sprintf("alpha1=%f", alpha1));

    # We require that the NB2 dispersion at mu=100 is less than 1
    # i.e., exp(alpha0) (100/N)^alpha1 < 1
#    alpha0.upper = -alpha1*log(100.0/N);

  #  if (print.level>1)
    #  message(sprintf("Maximize l(alpha0, alpha1) on (%f, %f)",
      #                alpha0.lower, alpha0.upper));

   # obj = optimize(l.alpha0, interval=c(alpha0.lower, alpha0.upper),
     # alpha1=alpha1, tol = tol.alpha0, maximum=TRUE);

   # if (print.level>1)
     # message(sprintf("alpha0=%f, l=%f", obj$maximum, obj$objective));

   # obj$objective
 # }

 # Maximize the pl of alpha1
 # if (is.na(alpha1)) {
 #   alpha1 = optimize(pl.alpha1, interval=alpha1.bounds,
 #     tol=tol.alpha1, maximum=TRUE)$maximum;
#  }

  #  if (print.level>0) message(sprintf("alpha1=%f", alpha1));

#  alpha0.upper = -alpha1*log(100.0/N);
#  if (print.level>1)
#    message(sprintf("Maximize l(alpha0, alpha1) on (%f, %f)",
#                    alpha0.lower, alpha0.upper));

#  obj = optimize(l.alpha0, interval=c(-20, alpha0.upper),
#    alpha1=alpha1, tol = tol.alpha0, maximum=TRUE);

#  alpha0=obj$maximum;
#  l = obj$objective;

#  if (print.level>0)
#    message(sprintf("alpha0=%f, alpha1=%f, l=%f", alpha0, alpha1, l));

#  m = dim(y)[1];
#  n = dim(y)[2];
#  mu = matrix(0, m, n);
  
#  for (i in 1:m) {
#    obj = apl.phi.alpha.1(exp(alpha0), alpha1, y[i,], lib.sizes, x, beta,
#      print.level=print.level-1);
##  l[i] = obj$likelihood;
#    mu[i,] = obj$mu;
#  }

#  p = mu / (matrix(1, m, 1) %*% matrix(lib.sizes, 1, n));
#  phi=exp(alpha0)*p^alpha1;
#  list(phi=phi, mu=mu, e=p, lib.sizes=lib.sizes, alpha0=alpha0, alpha1=alpha1, l=l);

  m <- dim(y)[1];
  n <- dim(y)[2];
  mu <- matrix(0,m,n);
  # Estimate phi using estimators of  alpha0, alpha1 from real data 
  alpha0 = obs.fit$alpha0.est
  alpha1 = obs.fit$alpha1.est
#   alpha0 <- fit$models[[1]]$alpha0;
#   alpha1 <- fit$models[[1]]$alpha1;
    
  for (i in 1:m) {
    obj <- apl.phi.alpha.1(exp(alpha0),alpha1, y[i,], lib.sizes=lib.sizes, x=x, beta=beta,
         print.level=print.level-1);
    l[i] <- obj$likelihood;
    mu[i,] <- obj$mu;
  }
  
  p <- mu / (matrix(1,m,1)%*% matrix(lib.sizes, 1, n));
  phi <- exp(alpha0)*p^alpha1;
  list(phi=phi, mu=mu, e=p,lib.sizes=lib.sizes, alpha0=alpha0, alpha1=alpha1, l=l);
}
############################################################################################
