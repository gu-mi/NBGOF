
#' @title Empirical Probability Plot
#'
#' @description A plot method for empirical probability plot based on a "gofm" 
#' object obtained from the \code{\link{nb_gof_m}} outputs
#' 
#' @param x an object of class "gofm" (from \link{nb_gof_m} output)
#' @param conf.env confidence level for the envelope
#' @param data.note a note on how data are simulated or the source of data
#' @param leg.cex size of the legend fonts
#' @param put.lab logical value: whether to put gene index labels on the plot
#' @param ... 
#'
#' @rdname plot
#' @method plot gofm
#' @export
#' 
#' @details This is a generic function for plotting goodness-of-fit test results.
#' 
#' @return A plot
#' 
#' @author Gu Mi, Yanming Di, Dan Schafer
#' 
plot <- function(x, ...) {
 UseMethod("plot", x)
}

plot.gofm <- function(x, ..., conf.env=0.95, data.note="", leg.cex=1, put.lab = TRUE){
  
  ## quantities from a "gofm" object:
  model.fit <- x$model.fit
  counts.dim = x$counts.dim
  xx = x$design.mat
  pv.P <- x$pear.pval
  pv.T <- x$new.pval
  stat.sim <- x$stat.sim
  stat.sim.T <- x$stat.sim.T
  o.res.sim <- x$o.res.sim
  res.typic <- x$typ.res.sim  # for x-axis plotting
  res.obs <- x$o.res0         # for y-axis plotting
  stat0 <- x$stat0
  stat0.T <- x$stat0.T
  sim <- x$sim
  n.pts <- length(res.typic)      # number of points
  
  # determine what quantiles to use for plotting:
  alpha = 1-conf.env
  q.upp <- conf.env + (1-conf.env)/2
  q.low <- (1-conf.env)/2
  quant <- apply(o.res.sim, 2, quantile, probs=c(q.low,q.upp))
  
  ## plotting options:
  min.x <- min(res.typic)
  max.x <- max(res.typic)
  min.y <- min(quant[1, ])
  max.y <- max(quant[2, ])
  
  # estiamte p.hat:
  out.ind = (res.obs < quant[1, ]) | (res.obs > quant[2, ])
  nout = sum(out.ind)
  p.hat = nout/n.pts
  
  ## plot epp:
  plot(res.typic[!out.ind], res.obs[!out.ind], xlim = c(min.x, max.x), 
       ylim = c(min.y, max.y), xlab =  "Residuals from Simulations",
       ylab =  "Original Residuals",
       main = paste("Testing ", model.fit, " Model Fit", "\n",
                    " (sim = ", sim, ", CI = ", conf.env*100, "%)", 
                    sep=""),...)
  points(res.typic[out.ind], res.obs[out.ind], col="red", pch=".", cex=2)
  abline(a=0,b=1,col="blue",lty="dashed")
  
  # quantiles (but not confidence bands -- just points)
  points(res.typic, quant[1, ], col="cyan2",cex=2, pch=".")
  points(res.typic, quant[2, ], col="cyan2",cex=2, pch=".")
  if (put.lab == TRUE){
    labels = seq(1,n.pts)
    toLabel = ceiling(labels[out.ind]/dim(xx)[1])
    xouts = res.typic[out.ind]
    youts = res.obs[out.ind]
    textxy(xouts, youts, toLabel, m=c(mean(xouts),mean(youts)))
  }
  
  legend("topleft",bty="n", legend=c(
    paste("Pear.T.pval = ", pv.P),
    paste("New.T.pval = ", pv.T),
    paste("Pear.Stat = ", round(stat0,2)),
    paste("Dim.Count = ", counts.dim),
    paste("#outside = ", nout, " (",round(p.hat*100,2),"%)", sep="")),
         cex=leg.cex)   
  legend("bottomright", bty="n", legend=paste("Data:", data.note)) 
}