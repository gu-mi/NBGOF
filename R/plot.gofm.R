
#' @title Plot Methods for "gofm" Objects (Histograms and Q-Q Uniform Plots of P-values)
#'
#' @description A plot method for evaluating the p-values based on a "gofm" 
#' object obtained from the \code{\link{nb_gof_m}} outputs.
#' 
#' @param x an object of class "gofm" from the \code{\link{nb_gof_m}} output
#' @param model the model name used to fit the data (used for histogram title)
#' @param type either a quantile-quantile uniform plot (\code{qq}), or a histogram (\code{hist})
#' @param logscale whether use logscale for the uniform Q-Q plot (default is TRUE) as recommended
#' in the \code{\link{gap}} package
#' @param ... for future use
#' 
#' @method plot gofm
#' @rdname plot.gofm
#' @export
#' 
#' @importFrom gap qqunif
#' 
#' @details This is a generic function for plotting goodness-of-fit test results. Users can
#' specify either a quantile-quantile (\code{qq}) uniform plot, or a histogram (\code{hist})
#' of the GOF p-values.
#' 
#' @usage plot(x, model=NULL, type="qq", logscale=TRUE, ...)
#' 
#' @return A quantile-quantile uniform plot or a histogram of the GOF p-values.
#' 
#' @seealso \code{\link{nb_gof_m}} for simulated data examples.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
plot.gofm = function(x, model=NULL, type="qq", logscale=TRUE, ...){
  
  stopifnot(type %in% c("qq", "hist"))
  
  ## quantities from a "gofm" object:
  model = x$model
  counts.dim = x$counts.dim
  sim.size = x$sim
  v.pvals = x$v.pvals
  
  ## plot qq-plot or histogram of m p-values
  if (type == "qq"){
    qqunif(v.pvals, logscale=logscale, main=paste("Q-Q Uniform Plot of", model))
  }
  else if (type == "hist"){
    hist(v.pvals, main=paste("Model:", model), xlab="P-values")
  }
  
}
