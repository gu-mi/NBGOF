
## code from Yanming on 2012/11/08 (concise version of estimate.disp.1.R)

## Estimate dispersion for each gene separately

library(numDeriv);

irls.nb = function(y, s, x, phi, beta0, ..., print.level=0) {
  m = dim(y)[1];
  n = dim(y)[2];
  p = dim(x)[2];
  
  phi = matrix(phi, m, n, byrow=TRUE);
  
  if (print.level > 0)
    print("Estimating NB regression coefficients using IRLS.");
  
  res=list(mu=matrix(NA, m, n), beta=matrix(NA, m, p), phi=matrix(NA, m, n),
           conv=logical(m), iter=numeric(m));
  
  if (print.level > 1) {
    ## Set up progress bar
    pb=txtProgressBar(style=3);
  }
  
  for (i in 1:m) {
    if (print.level > 1) {
      setTxtProgressBar(pb, i/m);
    }   
    res0 = irls.nb.1(y[i,], s, x, phi[i,], beta0,
                     ...,
                     print.level=print.level-1);
    res$mu[i,] =  res0$mu;
    res$beta[i,] = res0$beta;
    res$conv[i] = res0$conv;
    res$iter[i] = res0$iter;
  }
  
  if (print.level>1) close(pb);
  
  res;
}


## For fitting NB2 regression model, when glm.nb() has convergence issues

irls.nb.1 <- function(y, s, x, phi, beta0=rep(NA,p),
                      maxit=50, tol.mu=1e-3/length(y), print.level=0) {
  
  nobs = as.integer(dim(x)[1]);
  p = as.integer(dim(x)[2]);
  
  ## Indices to fixed and free components of beta
  id1 = (1:p)[is.na(beta0)];  # NA beta's
  id0 = (1:p)[!is.na(beta0)]; # non-NA beta
  q = length(id0);  # number of beta's (non-NA) to test (currently =1)
  nvars = p - q;
  
  ## Offset
  beta = beta0;
  offset = matrix(x[, id0], nobs, q) %*% beta[id0];
  
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
    
    z = eta - offset + (y - mu)/mu;
    w = drop(mu/sqrt(varmu));
    
    ## call Fortran code to perform weighted least square
    ## cf. glm.fit()
    epsilon = 1e-7;
    fit = .Fortran("dqrls",
                   qr = x[,id1] * w, n = nobs,
                   p = nvars, y = w * z, ny = 1L,
                   tol = epsilon,
                   coefficients = double(nvars),
                   residuals = double(nobs),
                   effects = double(nobs),
                   rank = integer(1L),
                   pivot = 1L:nvars,
                   qraux = double(nvars),
                   work = double(2 * nvars),
                   PACKAGE = "base");
    
    if (any(!is.finite(fit$coefficients))) {
      warning(gettextf("non-finite coefficients at iteration %d", iter),
              domain = NA);
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


log.likelihood.nb= function(kappa, mu, y) {
  ## p = mu/(mu+kappa);
  ## l0= sum(log(gamma(kappa+y)) - log(gamma(kappa)) - log(gamma(1+y))
  ##       + y * log(mu) + kappa * log(kappa) - (y+kappa) * log(mu+kappa));
  
  sum(dnbinom(y, kappa, mu=mu, log=TRUE));
}

pl.phi.1 = function(phi, y, lib.sizes, xx,
                    beta0=rep(NA, dim(xx)[2]),
                    log.dispersion=TRUE,
                    print.level=0) {
  
  if (log.dispersion) {
    phi = exp(phi);
  }
  fit = irls.nb.1(y, lib.sizes, xx, phi, beta0);
  log.likelihood.nb(1/phi, fit$mu, y);
  
}

mle.phi.1 = function(likelihood, interval, y, s, x, beta0, 
                     information=TRUE,
                     log.dispersion=TRUE,
                     print.level=1) {
  
  ## Maximize the profile likelihood of phi
  obj = optimize(likelihood, interval=interval,
                 y=y, lib.sizes=s, xx=x, beta0=beta0,
                 log.dispersion=log.dispersion,
                 print.level=print.level-1, 
                 maximum=TRUE);
  
  if (information) {
    h = -hessian(likelihood, obj$maximum,
                 y=y, lib.sizes=s, xx=x, beta0=beta0,
                 log.dispersion=log.dispersion); 
  } else {
    h = NA;
  }
  
  ## Return the MLE
  list(maximum=obj$maximum, information=h);
}

estimate.disp.1 = function(counts, lib.sizes, x, beta0, 
                           likelihood, interval = c(-20, 7),
                           log.dispersion=TRUE,
                           information=TRUE,
                           print.level=0) {
  n = dim(counts)[1];
  
  phi = numeric(n);
  i.hat = numeric(n);
  
  ## Set up progress bar
  #pb=txtProgressBar(style=3);
  
  ## Estimate the disperison parameter for each gene separately
  for (i in 1:n) {
    #setTxtProgressBar(pb, i/n);
    res =  mle.phi.1(likelihood, interval, counts[i,], lib.sizes, x, beta0,
                     log.dispersion=log.dispersion,
                     information=information);
    phi[i]= res$maximum;
    i.hat[i] = res$information;
  }
  #close(pb);
  
  list(estimate=phi, information=i.hat);
}


nb.regression.1 = function(y, s, x, beta, likelihood=pl.phi.1) {
  
  if (!is.matrix(y)) {
    dim(y) = c(1, length(y));
  }
  
  ## By default, we estimate the MLE of log(phi)
  log.phi.hat = estimate.disp.1(y, s, x, beta,
                                information=FALSE,
                                likelihood=likelihood)$estimate;
  phi.hat = exp(log.phi.hat);
  
  ## Fit NB regression
  fit = irls.nb(y, s, x, phi.hat, beta0=beta);
  
  fit$phi = matrix(phi.hat, dim(y)[1], dim(y)[2]); 
  
  ## If you need variance, uncomment the following line
  fit$v = fit$mu + fit$phi * fit$mu^2;
  
  fit
}