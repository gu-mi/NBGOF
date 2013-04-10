
# This file is to be sourced for GOF test -- glm.nbp.1()
# phi0 does NOT have to be specified beforehand (bug fixed)

# source('utils.R');  # for printf() function
# library(compiler)

# log.likelihood.nb() will be used in apl.phi.alpha.1()
# Negative Binomial (NB2) formulation, with kappa taking care of NBP
log.likelihood.nb <- function(kappa, mu, y) {
  sum(dnbinom(y, size=kappa, mu=mu, log=TRUE));
  # dnbinom(x, size, prob, mu, log=T)
}

# n.full.likelihood() will be used in optim.full.likelihood()
n.full.likelihood <- function(par, y, s, x) {
  alpha0 = par[1]; # 1st para: alpha0 --> knowing alpha0 --> knowing phi0
  phi0 = exp(alpha0);
  alpha1 = par[2]; # 2nd para: alpha1 --> knowing alpha1 --> knowing alpha
  beta = par[-(1:2)];  # par[alpha0, alpha1, beta]
  mu = s * exp(x%*%beta);
  phi = phi0 * (mu/s)^alpha1; # knowing phi0 and alpha1 --> knowing phi
  -sum(dnbinom(y, size=1/phi, mu=mu, log=TRUE));
  # NB2 pdf, but phi --> as a function of mu and alpha1 --> accounts for NB
  # result: -logl
}

optim.full.likelihood <- function(y, s, x, trace=0) {
  # trace=0: used in optim(); positive numbers display tracing info.
  
  fit = irls.nb.1(y, s, x, 0.1);
  # phi0=0.1 for initials of beta.hat's from NB2 modeling
  # phi = phi0 for NB2 case; also, var = mu + 0.1*mu^2
  # initial values --> pars[alpha0, alpha1, beta]
  # alpha0=log(0.1), since phi0=0.1; alpha1=0 --> NB2 case;
  # fit$beta --> fitted NB2 betas
  # data: y, s, x -- cf. n.full.likelihood(par,y,s,x)
  
  obj = optim(c(log(0.1), 0, fit$beta), n.full.likelihood, y=y, s=s, x=x,
              control=list(trace=trace), method="BFGS");
  beta = obj$par[-(1:2)];
  mu = s * exp(x%*%beta);
  list(alpha0 = obj$par[1], alpha1 = obj$par[2], phi0=exp(obj$par[1]),
       beta=beta, mu=mu);
}

##' Estimate the regression coefficients in an NB GLM model with known
##' dispersion parameters
##'
##' This function estimate <beta> using iterative reweighted least
##' squares (IRLS) algorithm, which is equivalent to Fisher scoring.
##' We used the glm.fit code as a template.
##'
##' @title Estimate the regression coefficients in an NB GLM model
##' @param y an n vector of counts
##' @param s a scalar or an n vector of effective library sizes
##' @param x a n by p design matrix
##' @param phi a scalar or an n-vector of dispersion parameters
##' @param beta0 a vector specifying known and unknown components of
##' the regression coefficients: non-NA components are hypothesized
##' values of beta, NA components are free components
##' @param maxit 
##' @param tol.mu convergence criteria
##' @param print.level 
##' @return a list of the following components:
##'  beta, a p-vector of estimated regression coefficients
##'  mu, an n-vector of estimated mean values
##'  converged, logical. Was the IRLS algorithm judged to have converged?
##'  @useDynLib NBGOF Cdqrls
##'  @keywords internal
irls.nb.1 = function(y, s, x, phi, beta0=rep(NA,p),
                     maxit=50, tol.mu=1e-3/length(y), print.level=0) {
  
  nobs = as.integer(dim(x)[1]);
  p = as.integer(dim(x)[2]);
  
  ## Indices to fixed and free components of beta
  id1 = (1:p)[is.na(beta0)];
  id0 = (1:p)[!is.na(beta0)];
  q = length(id0);
  nvars = p - q;
  
  
  ## Offset
  beta = beta0;
  offset = matrix(x[, id0], nobs, q) %*% beta[id0];
  
  ## Initial estiamte of beta
  ## eta = log((y+0.5)/s);
  ## beta[id1] = qr.solve(x[, id1], eta - offset);
  ## eta = drop(x %*% beta);
  mu = y + (y==0)/6;
  eta = log(mu/s);
  
  ##------------- THE Iteratively Reweighting L.S. iteration -----------
  conv = FALSE;
  for (iter in 1L:maxit) {
    varmu = mu + phi * mu^2;
    if (any(is.na(varmu)))
      stop("NAs in V(mu)")
    if (any(varmu == 0))
      stop("0s in V(mu)")
    
    ## mu.eta.val = mu;
    ## if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")
    ## drop observations for which w will be zero
    ## good = (mu.eta.val != 0)
    z = eta - offset + (y - mu)/mu;
    w = drop(mu/sqrt(varmu));
    
    ## call Fortran code to perform weighted least square
    epsilon = 1e-7;
    
    ## call Fortran code via C wrapper
    fit = .Call(Cdqrls, x[, id1, drop=FALSE] * w,  w * z, epsilon);
    
    ##    fit = .Fortran("dqrls",
    ##                    qr = x[,id1] * w, n = nobs,
    ##                    p = nvars, y = w * z, ny = 1L,
    ##                    tol = epsilon,
    ##                    coefficients = double(nvars),
    ##                    residuals = double(nobs),
    ##                    effects = double(nobs),
    ##                    rank = integer(1L),
    ##                    pivot = 1L:nvars,
    ##                    qraux = double(nvars),
    ##                    work = double(2 * nvars),
    ##                    PACKAGE = "base");
    
    if (any(!is.finite(fit$coefficients))) {
      warning(gettextf("non-finite coefficients at iteration %d", iter), domain = NA);
      break
    }
    
    ## stop if not enough parameters
    if (nobs < fit$rank)
      stop(gettextf("X matrix has rank %d, but only %d observations",
                    fit$rank, nobs), domain = NA)
    
    ## calculate updated values of eta and mu with the new coef:
    muold = mu;
    beta[id1[fit$pivot]] = fit$coefficients;
    eta = drop(x %*% beta);
    mu = s*exp(eta);
    
    ## check for convergence of mu
    ## ! when mu converges, beta may not
    if (max(abs(mu - muold)) < tol.mu) {
      conv = TRUE;
      break
    }
  }
  
  ##-------------- end IRLS iteration -------------------------------
  list(mu=mu, beta=beta, conv=conv, iter=iter);
}


# irls.nb = function(y, s, x, phi, beta0, ..., print.level=0) {
#   m = dim(y)[1];
#   n = dim(y)[2];
#   p = dim(x)[2];
#   
#   phi = matrix(phi, m, n, byrow=TRUE);
#   
#   if (print.level > 0)
#     print("Estimating NB regression coefficients using IRLS.");
#   
#   res=list(mu=matrix(NA, m, n), beta=matrix(NA, m, p),
#            conv=logical(m), iter=numeric(m));
#   
#   if (print.level > 1) {
#     ## Set up progress bar
#     pb=txtProgressBar(style=3);
#   }
#   
#   for (i in 1:m) {
#     if (print.level > 1) {
#       setTxtProgressBar(pb, i/m);
#     }
#     
#     res0 = irls.nb.1(y[i,], s, x, phi[i,], beta0,
#                      ...,
#                      print.level=print.level-1);
#     res$mu[i,] =  res0$mu;
#     res$beta[i,] = res0$beta;
#     res$conv[i] = res0$conv;
#     res$iter[i] = res0$iter;
#   }
#   
#   if (print.level>1) close(pb);
#   
#   res;
# }

##' Estimate the regression coefficients in an NBP GLM model for one gene
##'
##' This function estimate <beta> using iterative reweighted least
##' squares (IRLS) algorithm, which is equivalent to Fisher scoring.
##' We used the glm.fit code as a template.
##'
##' Note that we will igore the dependence of the dispersion parameter
##' (reciprical of the shape parameter) on beta. In other words, the
##' estimate is the solution to
##'
##'   dl/dmu dmu/dbeta = 0
##'
##' while we igored the contribution of 
##'
##'   dl/dkappa dkappa/dbeta
##'
##' to the score equation.
##'
##'
##' @title Estiamte the regression coefficients in an NBP GLM model
##' @param y an n vector of counts
##' @param s a scalar or an n vector of effective library sizes
##' @param x a n by p design matrix
##' @param phi0
##' @param alpha1  phi= phi0 (mu/s)^alpha1
##' @param beta0 the regression coefficients: non-NA components are hypothesized
##' values of beta, NA components are free components
##' @param tol.mu convergence criteria
##' @return a list of the following components:
##'  beta, a p-vector of estimated regression coefficients
##'  mu, an n-vector of estimated mean values
##'  converged, logical. Was the IRLS algorithm judged to have converged?
##'  @useDynLib NBGOF Cdqrls
##'  @keywords internal
irls.nbp.1 = function(y, s, x, phi0, alpha1, beta0=rep(NA, p),
                      maxit=50, tol.mu=1e-3/length(y), print.level=1) {
  
  nobs = as.integer(dim(x)[1]);
  p = as.integer(dim(x)[2]);
  
  ## Indices to fixed and free components of beta
  id1 = (1:p)[is.na(beta0)];
  id0 = (1:p)[!is.na(beta0)];
  q = length(id0);
  nvars = p - q;
  
  ## Offset
  beta = beta0;
  offset = matrix(x[, id0], nobs, q) %*% beta[id0];
  
  ## Initial estiamte of beta
  ## eta = log((y+0.5)/s);
  ## beta[id1] = qr.solve(x[, id1], eta - offset);
  ## eta = drop(x %*% beta);
  mu = y + (y==0)/6;
  eta = log(mu/s);
  
  ##------------- THE Iteratively Reweighting L.S. iteration -----------
  conv = FALSE;
  for (iter in 1L:maxit) {
    phi = phi0 * (mu/s)^alpha1;
    varmu = mu + phi * mu^2;
    if (any(is.na(varmu)))
      stop("NAs in V(mu)")
    if (any(varmu == 0))
      stop("0s in V(mu)")
    
    ## mu.eta.val = mu;
    ## if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")
    ## drop observations for which w will be zero
    ## good = (mu.eta.val != 0)
    z = eta - offset + (y - mu)/mu;
    w = drop(mu/sqrt(varmu));
    
    ## call Fortran code to perform weighted least square
    epsilon = 1e-7;
    
    ## call Fortran code via C wrapper
    fit = .Call(Cdqrls, x[, id1, drop=FALSE] * w,  w * z, epsilon);
    
    ##    fit = .Fortran("dqrls",
    ##                    qr = x[,id1] * w, n = nobs,
    ##                    p = nvars, y = w * z, ny = 1L,
    ##                    tol = epsilon,
    ##                    coefficients = double(nvars),
    ##                    residuals = double(nobs),
    ##                    effects = double(nobs),
    ##                    rank = integer(1L),
    ##                    pivot = 1L:nvars,
    ##                    qraux = double(nvars),
    ##                    work = double(2 * nvars),
    ##                    PACKAGE = "base");
    
    if (any(!is.finite(fit$coefficients))) {
      warning(gettextf("non-finite coefficients at iteration %d", iter), domain = NA);
      break
    }
    
    ## stop if not enough parameters
    if (nobs < fit$rank)
      stop(gettextf("X matrix has rank %d, but only %d observations",
                    fit$rank, nobs), domain = NA)
    
    ## calculate updated values of eta and mu with the new coef:
    muold = mu;
    beta[id1[fit$pivot]] = fit$coefficients;
    eta = drop(x %*% beta);
    mu = s*exp(eta);
    
    ## check for convergence of mu
    ## ! when mu converges, beta may not
    if (max(abs(mu - muold)) < tol.mu) {
      conv = TRUE;
      break
    }
  }
  
  ##-------------- end IRLS iteration -------------------------------
  
  phi = phi0 * (mu/s)^alpha1;
  list(mu=mu, beta=beta, phi=phi, conv=conv, iter=iter);
}


# irls.nbp = function(y, s, x, phi0, alpha1,
#                     maxit=50, tol.mu=0.01, print.level=0) {
#   
#   n = dim(y)[1];
#   K = dim(x)[2];
#   
#   if (print.level > 0)
#     print("Estimating NB regression coefficients using IRLS.");
#   
#   res=list(mu=y, beta=matrix(0, n, K), conv=rep(FALSE, n), iter=numeric(n));
#   
#   if (print.level > 0) pb = txtProgressBar(style=3);
#   
#   for (i in 1:n) {
#     if (print.level > 0) {
#       setTxtProgressBar(pb, i/n);
#     }
#     res0 = irls.nbp.1(y[i,], s, x, phi0, alpha1, maxit=maxit,
#                       tol.mu=tol.mu, print.level=print.level-1);
#     res$mu[i,] =  res0$mu;
#     res$beta[i,] = res0$beta;
#     res$conv[i] = res0$conv;
#     res$iter[i] = res0$iter;
#   }
#   
#   if (print.level>0) close(pb);
#   
#   res;
# }


## -------------------------------------------------------------------------
## Log adjusted profile likelihood (APL) of parameters in dispersion model:
## phi = phi0*p^alpha1
## -------------------------------------------------------------------------
apl.phi.alpha.1 <- function(phi0, alpha1, y, lib.sizes, x, beta0,
                            print.level=0) {
  
  ## Estimate regresssion coefficients in the NBP regression model
  fit = irls.nbp.1(y, lib.sizes, x, phi0, alpha1, beta0);
  
  kappa = 1/fit$phi;  # "size" in log.likelihood.nb()
  
  ## Compute the information matrices
  mu.hat = fit$mu;
  v.hat = drop(mu.hat + fit$phi * mu.hat^2);
  j.hat <- t(x) %*%
    diag(mu.hat^2 * (y/mu.hat^2 - (y+kappa)/(mu.hat+kappa)^2) -
           (y-mu.hat)*mu.hat/v.hat) %*% x;
  
  ## The likelihood of (kappa, mu) (different from the l.hat in
  ## irls.nb.1, which is the likelihood of mu).
  l.hat = log.likelihood.nb(kappa, fit$mu, y);
  
  if (print.level>2) {
    print(list(mu.hat=mu.hat, v.hat=v.hat, l.hat=l.hat, j.hat=j.hat));
  }
  l.hat - 0.5 * log(det(j.hat));
  # equation (8) in paper -- adjusted profile likelihood
}

## -------------------------------------------------------------------------
## Estimate the NBP dispsersioon model and fit NBP regression model
## -------------------------------------------------------------------------
glm.nbp.1 <- function(y, s, x,
                      beta0 = rep(NA, dim(x)[2]),
                      alpha1=NA,
                      mu.lower=1,
                      mu.upper=Inf,
                      tol.alpha0 = 1e-4,
                      tol.alpha1 = 1e-3,
                      print.level=3) {
  
  ## Mean relative counts
  n = length(y);
  mu.lower = max(mu.lower, 1);
  
  if (print.level>0) {
    message("Estimating the dispersion model:");
    message("log(dispersion) = alpha0 + alpha1 log(mu/lib.size)");
  }
  
  ## Bounds for alpha0 and alpha1.
  ## We are mainly interested in alpha \in [0.5, 2.5]
  ## --> alpha1 \in [-1.5, 0.5]
  
  ## alpha.bounds = c(1 - tol.alpha, 2 + tol.alpha);  # tol.alpha = 0.5
  alpha1.bounds = c(-1.5, 0.5);
  # 2 + alpha1 = alpha, so alpha \in [0.5, 2.5]
  alpha0.lower = -20;
  # alpha0 = log(phi0)  ---> phi0.lower = 2e-9 (+ over-disp)
  
  N = max(s);
  
  # l.alpha0() --> to be passed to optimize()
  l.alpha0 <- function(alpha0, alpha1) {
    if (print.level>2) {
      printf("alpha0=%f, phi=%f when mu=1000 in a library of size %.0f.\n",
             alpha0, exp(alpha0)*(1000/N)^alpha1, N);
    }
    apl.phi.alpha.1(exp(alpha0), alpha1, y, s, x, beta0,
                    print.level=print.level-1);
  }
  
  ## The adjusted profile likelihood of alpha
  pl.alpha1 <- function(alpha1) {
    if (print.level>0) message(sprintf("alpha1=%f", alpha1));
    
    ## We require that the NB2 dispersion at mu=100 is less than 1
    ## ow, var = mu + phi*mu^2 ~= 10000*phi, std > 100 = mu -- unrealistic!
    ## i.e., phi = exp(alpha0)*(100/N)^alpha1 < 1, thus...
    alpha0.upper = -alpha1*log(100.0/N);
    
    if (print.level>1)
      message(sprintf("Maximize l(alpha0, alpha1) on (%f, %f)",
                      alpha0.lower, alpha0.upper));
    
    # now maximize l.alpha0 for fixed alpha1
    obj = optimize(l.alpha0, interval=c(alpha0.lower, alpha0.upper),
                   alpha1=alpha1, tol = tol.alpha0, maximum=TRUE);
    
    if (print.level>1)
      message(sprintf("alpha0=%f, l=%f", obj$maximum, obj$objective));
    
    obj$objective
  }
  
  ## Maximize the pl of alpha1
  if (is.na(alpha1)) {
    alpha1 = optimize(pl.alpha1, interval=alpha1.bounds,
                      tol=tol.alpha1, maximum=TRUE)$maximum;
  }
  if (print.level>0) message(sprintf("alpha1=%f", alpha1));
  
  alpha0.upper = -alpha1*log(100.0/N);
  if (print.level>1)
    message(sprintf("Maximize l(alpha0, alpha1) on (%f, %f)",
                    alpha0.lower, alpha0.upper));
  
  obj = optimize(l.alpha0, interval=c(alpha0.lower, alpha0.upper),
                 alpha1=alpha1, tol = tol.alpha0, maximum=TRUE);
  alpha0=obj$maximum;
  l = obj$objective;
  
  if (print.level>0)
    message(sprintf("alpha0=%f, alpha1=%f, l=%f", alpha0, alpha1, l));
  
  fit = irls.nbp.1(y, s, x, exp(alpha0), alpha1, beta0);
  
  fit$alpha0 = alpha0;
  fit$alpha1 = alpha1;
  p = fit$mu/s;
  fit$phi=exp(alpha0)*p^alpha1;
  fit$l = l;
  # Pearson's residuals:
  fit$p.res = (y - fit$mu) / sqrt( fit$mu + fit$phi * fit$mu^2 );
  # Pearson's GOF test statistic:
  # fit$p.stat = sum((fit$p.res)^2);
  
  fit;
}