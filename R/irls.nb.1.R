
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
