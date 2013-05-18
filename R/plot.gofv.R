
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
#' @author Gu Mi, Yanming Di, Daniel Schafer
#'  
plot.gofv = function(x, 
                     conf.env=0.95, 
                     data.note=NULL, 
                     leg.cex=1, 
                     ...){
  
  ## quantities from a "gofv" object:
  model.fit = x$model.fit
  samp.size = x$samp.size
  pv.D = x$new.pval
  pv.P = x$pear.pval
  sim = x$sim
  
  ord.res.sim.mat = x$ord.res.sim.mat
  ord.typ.res.sim = x$ord.typ.res.sim         # for x-axis plotting
  ord.res.obs = x$ord.res.sim.mat[(sim+1), ]  # for y-axis plotting
  
  n.pts = length(ord.res.obs)  # number of points (samples)
  
  # determine what quantiles to use for plotting:
  alpha = 1-conf.env
  q.upp = conf.env + alpha/2
  q.low = alpha/2
  quant = apply(ord.res.sim.mat[-(sim+1), ], 2, quantile, probs = c(q.low, q.upp))
  
  ## plotting options:
  min.x = min(ord.typ.res.sim)
  max.x = max(ord.typ.res.sim)
  min.y = min(quant[1, ], ord.res.obs)
  max.y = max(quant[2, ], ord.res.obs)
  
  # estiamte p.hat:
  out.ind = ( (ord.res.obs < quant[1, ]) | (ord.res.obs > quant[2, ]) )
  nout = sum(out.ind)
  p.hat = nout/n.pts
  
  ## plot epp:
  if (model.fit == "NB2"){
    plot(ord.typ.res.sim[!out.ind], ord.res.obs[!out.ind], 
         xlim = c(min.x, max.x), 
         ylim = c(min.y, max.y), 
         xlab =  "Medians of Ordered MC Residuals (NB2)",
         ylab =  "Ordered Pearson Residuals",
         main = paste("Testing", model.fit, "Regression Fit"),
         ...)
  }
  else if (model.fit == "NBP"){
    plot(ord.typ.res.sim[!out.ind], ord.res.obs[!out.ind], 
         xlim = c(min.x, max.x), 
         ylim = c(min.y, max.y), 
         xlab =  "Medians of Ordered MC Residuals (NBP)",
         ylab =  "Ordered Pearson Residuals",
         main = paste("Testing", model.fit, "Regression Fit"),
         ...)
  }
  
  
  # highlight points outside:
  points(ord.typ.res.sim[out.ind], ord.res.obs[out.ind], col="red", ...)
  abline(a=0,b=1,col="blue",lty="dashed")
  
  # quantiles for epp:
  points(ord.typ.res.sim, quant[1, ], type="l", col="blue")
  points(ord.typ.res.sim, quant[2, ], type="l", col="blue")
  
  if (pv.D > 0.01 & pv.P > 0.01){
    legend("topleft",bty="n", legend=c(
      #paste("# sample = ", samp.size),
      #paste("# pts. outside envelope = ", nout, " (",round(p.hat*100,2),"%)", sep=""),
      paste("MC GOF p-value (Pear.Stat.) =", round(pv.P,2)),
      paste("MC GOF p-value (Sq.Vert.Dist.) =", round(pv.D,2))),
           cex=leg.cex)
  }
  #
  if (pv.D < 0.005 & pv.P > 0.01){
    legend("topleft",bty="n", legend=c(
      #paste("# sample = ", samp.size),
      #paste("# pts. outside envelope = ", nout, " (",round(p.hat*100,2),"%)", sep=""),
      paste("MC GOF p-value (Pear.Stat.) =", round(pv.P,2)),
      paste("MC GOF p-value (Sq.Vert.Dist.) < 0.01")),
           cex=leg.cex)
  }
  #
  if (pv.D > 0.01 & pv.P < 0.005){
    legend("topleft",bty="n", legend=c(
      #paste("# sample = ", samp.size),
      #paste("# pts. outside envelope = ", nout, " (",round(p.hat*100,2),"%)", sep=""),
      paste("MC GOF p-value (Pear.Stat.) < 0.01"),
      paste("MC GOF p-value (Sq.Vert.Dist.) =", round(pv.D,2))),
           cex=leg.cex)
  } 
  #
  if (pv.D < 0.005 & pv.P < 0.005){
    legend("topleft",bty="n", legend=c(
      #paste("# sample = ", samp.size),
      #paste("# pts. outside envelope = ", nout, " (",round(p.hat*100,2),"%)", sep=""),
      paste("MC GOF p-value (Pear.Stat.) < 0.01"),
      paste("MC GOF p-value (Sq.Vert.Dist.) < 0.01")),
           cex=leg.cex)
  }  
  ##  
#   legend("bottomright", bty="n", legend=c(data.note,
#     paste("(sim size = ", sim, ", CI = ", conf.env*100, "%)", sep="")),
#          cex=leg.cex)
  legend("bottomright", bty="n", legend=data.note, cex=leg.cex)
}