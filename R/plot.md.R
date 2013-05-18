

#' @title mean-dispersion plot (log-log) scale for paper use
#' 
#' @description This function is designed to 
#' 
#' @param y 
#' @param x an n-by-p design matrix.
#' @param model a string of characters specifying the negative binomial model used 
#' to fit the data. Currently the following dispersion models are available to be checked
#' for goodness-of-fit
#' @param scatter 
#' @param legend
#' @param data.type
#' 
#' @return An object
#' 
#' @details When the response is a vector, we can use this function to test
#' the goodness-of-fit of a specified negative binomial regression model. 
#' 
#' @usage 
#' plot.md(y, x, model = NULL, scatter = FALSE, legend = FALSE, data.type = NULL, ...)
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
#' 
plot.md = function(y, x, model = NULL, scatter = FALSE, legend = FALSE, data.type = NULL, ...){
  
  phi.line = function (mu, v, alpha = 2, ...) 
  {
    id = order(mu)
    mu = mu[id]
    v = v[id]
    phi = (v - mu)/mu^alpha
    #id = (phi > 0) & (mu > 0)
    #mu = mu[id]
    #phi = phi[id]
    if (length(mu) > 1000) {
      id = c(seq(1, length(mu) - 1, length = 1000), length(mu))
    }
    lines(mu[id], phi[id], ...)
    invisible()
  }
  
  # mean and variance and phi estimated from data:
  v = apply(y, 1, var)
  mu = apply(y, 1, mean)
  phi = (v - mu)/mu^2
  id = (phi > 0) & (mu > 0)
  total.pts = length(v)
  on.plot = sum(id)
  mu = mu[id]
  phi = phi[id]
  
  ## begin plotting 
  if (scatter){
    plot(mu, phi, log="xy", pch=".", cex=4, 
         xlab="Average Number of Reads (log)", ylab="Estimated NB Dispersion Parameter (log)",
         ...)
  }
  
  else if (!scatter){
    
    if (model == "NBP"){
      nb.data1 = prepare.nb.data(counts=y)
      nbpf = estimate.dispersion(nb.data1, x, print.level = 0)
      lines(nbpf$models[[1]]$mu[,1], nbpf$models[[1]]$phi[,1], col="red", lwd=2)
    }
    
    else if (model == "common"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      #
      e.com = estimateGLMCommonDisp(d, design, verbose = FALSE)
      com.fit = glmFit(d, design, dispersion = e.com$common.dispersion)
      mu.edgeR = com.fit$fitted.values
      phi.hat.m = com.fit$dispersion
      v.edgeR = mu.edgeR + phi.hat.m * mu.edgeR^2
      line.edgeR = list(mu = mu.edgeR, v = v.edgeR)
      #
      phi.line(mu.edgeR, v.edgeR, col="blue", alpha=2, lwd=2, lty=2)      
    }
    
    else if (model == "trended"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      #
      e.trd = estimateGLMTrendedDisp(d, design, min.n=100)
      trd.fit = glmFit(d, design, dispersion = e.trd$trended.dispersion)
      mu.edgeR = trd.fit$fitted.values
      phi.hat.m = trd.fit$dispersion
      v.edgeR = mu.edgeR + phi.hat.m * mu.edgeR^2
      #
      phi.line(mu.edgeR, v.edgeR, col="cyan", alpha=2, lwd=2, lty=5)
    }
    
    else if (model == "tagwise"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      #
      e.com = estimateGLMCommonDisp(d, design, verbose = FALSE)
      e.tag = estimateGLMTagwiseDisp(e.com, design) 
      #e.trd = estimateGLMTrendedDisp(d, design, min.n=25)
      #e.tag = estimateGLMTagwiseDisp(e.trd, design)
      tag.fit = glmFit(d, design, dispersion = e.tag$tagwise.dispersion)
      mu.edgeR = tag.fit$fitted.values
      phi.hat.m = tag.fit$dispersion
      v.edgeR = mu.edgeR + phi.hat.m * mu.edgeR^2
      #
      phi.line(mu.edgeR, v.edgeR, col="magenta", alpha=2, lwd=2, lty=4);
    }
  }
  
  if (legend){
    legend("bottomleft", bty="n", legend=c("common", "NBP", "trended", "tagwise"),
           lty=c(2,1,5,4), col=c("blue","red","cyan","magenta"),
           lwd=rep(2,4))
#     legend("topright", bty="n", legend=c(paste("plot",on.plot,"out of",total.pts,"points", 
#                                                sep=" "),
#                                          data.type = data.type))
  }
}

