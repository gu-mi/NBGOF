

#' @title The Mean-Dispersion Plot (on Log-log Scale) with Fitted Curves
#' 
#' @description This function provides a quick-and-dirty mean-dispersion plot (log-log scale)
#' with relative mean frequencies on the x-axis and estimated NB dispersions on the y-axis.
#' Several fitted curves from NB dispersion models are superimposed on the plot.
#' 
#' @param y an m-by-n count matrix of non-negative integers. For a typical
#' RNA-Seq experiment, this is the read counts with m genes and n samples.
#' @param x an n-by-p design matrix.
#' @param model a string of characters specifying the negative binomial model used 
#' to fit the data. Currently the following dispersion models are implemented: common, NBP,
#' trended (non-parametric), tagwise-common and tagwise-trended.
#' @param scatter logical: plot the points (\code{TRUE} if first invoking the plot). If 
#' \code{TRUE}, then a scatter plot is drawn, and the fitted NBP curve is superimposed.
#' @param legend logical: put the legend of dispersion models on plot (\code{TRUE} when plotting
#' the last curve).
#' @param ... for future use.
#' 
#' @return A mean-dispersion plot (log-log scale) with fitted curves from the following NB
#' dispersion models: 
#' common, NBP, trended (non-parametric), tagwise-common and tagwise-trended.
#' 
#' @details \strong{This function was originally used solely for paper figure 
#' productions. In the final release of the R package, it may be made invisible to end-users.
#' Also, current implementations are rather inflexible.}
#' 
#' @usage 
#' plot.md(y, x, model = NULL, scatter = FALSE, legend = FALSE, ...)
#' 
#' @seealso The Examples section of the \code{\link{arab}} dataset.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @export
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
plot.md = function(counts, x, model = NULL, scatter = FALSE, legend = FALSE, ...){
  
  ## begin plotting the points and NBP fitted curve
  if (scatter){
    m = dim(counts)[1]
    n = dim(counts)[2]
    nb.data = prepare.nb.data(counts = counts, 
                              lib.sizes = colSums(counts), 
                              norm.factors = rep(1, n))  
    nb.disp = estimate.dispersion(nb.data = nb.data, 
                                  x = x, 
                                  model = "NBP",
                                  method = "MAPL")    
    nbreg.fit = NBPSeq:::estimate.disp.mapl.nbp(y=counts, lib.sizes=colSums(counts), x=x)
    mu.hat = nbreg.fit$mu
    phi.hat = (1.5*rowSums((counts - mu.hat)^2) - rowSums(mu.hat))/rowSums(mu.hat^2)
    re.freq = (mu.hat / (matrix(1, m, 1) %*% matrix(colSums(counts), 1, n)))[ ,1]   
    id = (phi.hat > 0)  # may discard some phi.hat here not in the plotting
    #
    # scatterplot
    #
    plot(re.freq[id], phi.hat[id], log="xy", pch=".",
         xlab="Estimated Relative Frequency (log)", 
         ylab="Estimated NB Dispersion Parameter (log)",
         ...)
    lines(re.freq[id], nbreg.fit$phi[id,1], col="red", lwd=2, lty=1) 
  }
  #
  else if (!scatter){   
    
    #t(t(com.fit$fitted.values)/colSums(y))[,1]   # edgeR's pi
    
    if (model == "edgeR-common"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      #
      e.com = estimateGLMCommonDisp(d, design, verbose = FALSE)
      com.fit = glmFit(d, design, dispersion = e.com$common.dispersion)
      abline(h = com.fit$dispersion, lty=2, col="blue", lwd=2)
    }
    #
    else if (model == "edgeR-trended"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      
      ## debug(estimateGLMTrendedDisp);
      ## undebug(estimateGLMTrendedDisp);
      
      e.trd = estimateGLMTrendedDisp(d, design, min.n=100)
      ## trd.fit = glmFit(d, design, dispersion = e.trd$trended.dispersion)      
      trd.fit = glmFit(d, design, dispersion = 0.05);
      
      ## Gu, by calling estimateGLMTrendedDisp with a matrix (rather
      ## than a DGEList), you will get the estimated abundance
      ## offset = getOffset(d);
      ##
      ## disp = estimateGLMTrendedDisp(y, design, min.n=100, offset = offset);
      ## cor(disp$dispersion, e.trd$trended.dispersion);
      ## plot(disp$abundance,  disp$dispersion);
      ## lines(disp$abundance/1e6, disp$dispersion);
      
      ## YD
      id = order(trd.fit$fitted.values[,1]);
      pi.hat = trd.fit$fitted.values[id,1]/sum(y[,1]);
      phi.hat = e.trd$trended.dispersion[id];
      lines(pi.hat, phi.hat, col="cyan", lwd=3, lty=6)
      
      ## lines(sort(t(t(trd.fit$fitted.values)/colSums(y))[,1]), 
      ##      trd.fit$dispersion[order(t(t(trd.fit$fitted.values)/colSums(y))[,1])],
      ##      col="cyan", lwd=3, lty=6)
    }
    #
    else if (model == "edgeR-tagcom"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      #
      e.com = estimateGLMCommonDisp(d, design, verbose = FALSE)
      e.tgc = estimateGLMTagwiseDisp(e.com, design) 
      tgc.fit = glmFit(d, design, dispersion = e.tgc$tagwise.dispersion)
      lines(sort(t(t(tgc.fit$fitted.values)/colSums(y))[,1]), 
            tgc.fit$dispersion[order(t(t(tgc.fit$fitted.values)/colSums(y))[,1])],
            col="magenta", lwd=2, lty=4)
    }
    #
    else if (model == "edgeR-tagtrd"){
      grp.ids = factor(apply(x, 1, function(x) {
        paste(rev(x), collapse = ".")
      }), labels = seq(ncol(x)))
      d = DGEList(counts = y, lib.size = colSums(y), group = grp.ids)
      design = matrix(as.numeric(as.character(grp.ids)))
      #
      e.trd = estimateGLMTrendedDisp(d, design, verbose = FALSE)
      e.tgt = estimateGLMTagwiseDisp(e.trd, design) 
      tgt.fit = glmFit(d, design, dispersion = e.tgt$tagwise.dispersion)
      lines(sort(t(t(tgt.fit$fitted.values)/colSums(y))[,1]), 
            tgt.fit$dispersion[order(t(t(tgt.fit$fitted.values)/colSums(y))[,1])],
            col="yellow", lwd=2, lty=5)
    }
  }
  #
  if (legend){
    legend("bottomleft", bty="n", 
           legend=c("common", "NBP", "trended", "tagwise common", "tagwise trended"),
           lty=c(2,1,6,4,5), col=c("blue","red","cyan","magenta","yellow"),
           lwd=c(2,2,3,2,2))
  }
}