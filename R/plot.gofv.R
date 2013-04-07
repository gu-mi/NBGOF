
#' @title Plot Methods for "gofv" Objects
#'
#' @description A plot method for empirical probability plot based on a "gofv" 
#' object obtained from the \code{\link{nb_gof_v}} outputs
#' 
#' @param x an object of class "gofv" (from \link{nb_gof_v} output)
#' @param conf.env confidence level for the envelope
#' @param data.note a note on how data are simulated or the source of data
#' @param leg.cex size of the legend fonts
#' @param ... for future use
#' @method plot gofv
#' @rdname plot.gofv
#' @export
#' 
#' @details This is a generic function for plotting goodness-of-fit test results.
#' @usage plot(x, conf.env=0.95, data.note=NULL, leg.cex=1, ...)
#' 
#' @return An empirical probability plot with points outside a pre-specified confidence
#' envelope colored in red.
#' 
#' @seealso See \code{\link{nb_gof_v}} for simulated data examples
#' 
#' @author Gu Mi, Yanming Di, Dan Schafer
#'  
plot.gofv <- function(x, conf.env=0.95, data.note=NULL, leg.cex=1, ...){
  
  ## quantities from a "gofv" object:
  model.fit <- x$model.fit
  samp.size = x$samp.size
  pv.P <- x$pear.pval
  #pv.T <- x$new.pval
  #stat.sim <- x$stat.sim
  #stat.sim.T <- x$stat.sim.T
  o.res.sim <- x$o.res.sim
  res.typic <- x$typ.res.sim  # for x-axis plotting
  res.obs <- x$o.res0         # for y-axis plotting
  #stat0 <- x$stat0
  #stat0.T <- x$stat0.T
  sim <- x$sim
  n.pts <- length(res.typic)  # number of points
  
  # determine what quantiles to use for plotting:
  q.upp <- conf.env + (1-conf.env)/2
  q.low <- (1-conf.env)/2
  quant <- apply(o.res.sim, 2, quantile, probs=c(q.low,q.upp))
  
  ## plotting options:
  min.x <- min(res.typic)
  max.x <- max(res.typic)
  min.y <- min(quant[1, ])
  max.y <- max(quant[2, ])
  
  ## plot epp:
  out.ind = (res.obs < quant[1, ]) | (res.obs > quant[2, ])
  nout = sum(out.ind)
  plot(res.typic[!out.ind], res.obs[!out.ind], xlim = c(min.x, max.x), 
       ylim = c(min.y, max.y), xlab =  "Expected Ordered Res. from NB (via Simulation)",
       ylab =  "Ordered Pearson Residuals",
       main = paste("Testing ", model.fit, " Model Fit", "\n",
                    " (sim size = ", sim, ", CI = ", conf.env*100, "%)", 
                    sep=""),...)
  points(res.typic[out.ind], res.obs[out.ind], col="red", cex=5, pch=".")
  abline(a=0,b=1,col="blue",lty="dashed")
  
  # quantiles for epp:
  points(res.typic,quant[1, ], type="l")
  points(res.typic,quant[2, ], type="l")
  
  legend("topleft",bty="n", legend=c(
    paste("P-value (par.boot.Pear-stat) = ", pv.P),
    #paste("New.T.pval = ", pv.T),
    #paste("Pear.Stat = ", round(stat0,2)),
    paste("#sample = ", samp.size),
    paste("#outside envelope = ", nout, " (",round(nout/n.pts*100,2),"%)", sep="")),
         cex=leg.cex)
  legend("bottomright", bty="n", legend=paste("Actual Dist.:", data.note)) 
}