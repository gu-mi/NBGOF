
## For NB2 model fitting on the original & simulated datasets
# dispersion estimation: trended dispersion in edgeR

# update: 2012/11/30

model_edgeR_trended <- function(counts, x, lib.sizes = colSums(counts)){
  
#   nr = dim(counts)[1]
#   nc = dim(counts)[2]
#   n = nr * nc
  
  # initializations
#   mu.hat.m = matrix(0, nr = nr, nc = nc)
#   phi.hat.m = matrix(0, nr = nr, nc = nc)
#   res.m = matrix(0, nr = nr, nc = nc)
  
  ## edgeR tagwise dispersion:
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  design = model.matrix(~grp.ids, data=d$samples)
  e.com = estimateGLMCommonDisp(d, design, verbose=TRUE)
  e.trd = estimateGLMTrendedDisp(e.com, design)
  trd.fit = glmFit(d, design, dispersion=e.trd$trended.dispersion)
  
  # extract quantities:
  mu.hat.m = trd.fit$fitted.values   # mu may be close to 0
  phi.hat.m = trd.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
  res.m = (counts - mu.hat.m) / sqrt(v) 
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_trd_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
                         #fit.trd = trd.fit
  )
  return(model_trd_m_obj)
}