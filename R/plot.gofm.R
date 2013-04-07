
#' @title Plot Methods for "gofm" Objects
#'
#' @description A plot method for empirical probability plot based on a "gofm" 
#' object obtained from the \code{\link{nb_gof_m}} outputs
#' 
#' @param x an object of class "gofm" (from \link{nb_gof_m} output)
#' @param conf.env confidence level for the envelope
#' @param data.note a note on how data are simulated or the source of data
#' @param leg.cex size of the legend fonts
#' @param put.lab logical value: whether to put gene index labels on the plot
#' @param ... for future use
#' 
#' @method plot gofm
#' @rdname plot.gofm
#' @export
#' 
#' @importFrom calibrate textxy
#' 
#' @details This is a generic function for plotting goodness-of-fit test results.
#' @usage plot(x, conf.env=0.95, data.note=NULL, leg.cex=1, put.lab = TRUE, ...)
#' 
#' @return A diagnostic plot for evaluating model adequacy. Each gene index (row index) 
#' that falls outside a pre-specified confidence envelope is in red color and labels (if
#' \code{put.lab=TRUE}).
#' 
#' @seealso See \code{\link{nb_gof_m}} for simulated data examples
#' 
#' @author Gu Mi, Yanming Di, Dan Schafer
#' 
plot.gofm <- function(x, conf.env=0.95, data.note=NULL, leg.cex=1, put.lab = TRUE, ...){
  
  ## quantities from a "gofm" object:
  model.fit <- x$model.fit
  counts.dim = x$counts.dim
  xx = x$design.mat
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
       ylim = c(min.y, max.y), xlab =  "Expected Ordered Res. from NB (via Simulation)",
       ylab =  "Ordered Pearson Residuals",
       main = paste("Testing ", model.fit, " Model Fit", "\n",
                    " (sim size = ", sim, ", CI = ", conf.env*100, "%)", 
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
    paste("P-value (par.boot.Pear-stat) = ", pv.P),
    #paste("New.T.pval = ", pv.T),
    #paste("Pear.Stat = ", round(stat0,2)),
    paste("Count matrix dimension = ", counts.dim),
    paste("#outside envelope = ", nout, " (",round(p.hat*100,2),"%)", sep="")),
         cex=leg.cex)   
  legend("bottomright", bty="n", legend=paste("Actual Dist.:", data.note)) 
  
  invisible()
}


# #' Generic Plotting Function
# #' 
# #' @param model object of type gofv or gofm
# #' @param ... for future use
# #' @export 
# #' @return a plot
# #' 
# plot <- function(x, ...) {
#   UseMethod("plot", x)
# }

# #' Default Plotting Function
# #' 
# #' The default method for plot will return an error.  
# #' Since currently we have two types of response, vector and matrix, with classes
# #' "gofv" and "gofm" designed for effective plotting, it is not possible to write 
# #' a generic function to cater for any type of object. 
# #' 
# #' @param x an object
# #' @param ... for future use
# #' @method plot default
# #' @export 
# plot.default <- function(x, ...){
#   xx <- class(x)
#   stop(paste("No plot method defined for class", xx))
#   return(NULL)
# }
