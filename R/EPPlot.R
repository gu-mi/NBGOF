

#' @title Prepare Quantities for Empirical Probability Plots Using ggplot2
#' 
#' @description This function prepares quantities for making empirical probability plots (EPP) of negative binomial regression models. 
#' The result (a list) will be passed to the \code{\link{EPPlot}} function to make EPP.
#' 
#' @param gofv.obj an object of class "gofv" from the \code{\link{nb.gof.v}} output.
#' @param envelope the level of prediction band when drawing EPP.
#' @return a list of four components.
#' @keywords internal
#' 
eppdata = function(gofv.obj, envelope=0.95){
  
  if ( class(gofv.obj) != "gofv") {
    stop("You must pass an object of class 'gofv' as the first argument!")
  }
  
  model = gofv.obj$model
  sim = gofv.obj$sim
  pv.D = gofv.obj$new.pval
  pv.P = gofv.obj$pear.pval.1
  
  ord.res.sim.mat = gofv.obj$ord.res.sim.mat
  ord.typ.res.sim = gofv.obj$ord.typ.res.sim         # for x-axis plotting
  ord.res.obs = gofv.obj$ord.res.sim.mat[(sim+1), ]  # for y-axis plotting
  
  n.pts = length(ord.res.obs)  # number of points (samples)
  
  ## ------------------------------------------------------------------------------------------
  ## Simultaneous prediction interval
  ## ------------------------------------------------------------------------------------------
  col.ranks.mat = apply(gofv.obj$ord.res.sim.mat, 2, order)
  sv = numeric(sim)
  for (i in 1:sim){
    sv[i] = max(  max(col.ranks.mat[i,]), sim + 1 - min(col.ranks.mat[i,])  )
  }
  n.star = quantile(sv, probs=envelope)   # n* in Buja and Rolke, page 32, which determines the lower/upper simultaneous bounds
  ord.res.sim.mat2 = apply(gofv.obj$ord.res.sim.mat, 2, sort)
  sci.lower = ord.res.sim.mat2[sim + 1 - n.star, ]
  sci.upper = ord.res.sim.mat2[n.star, ]
  
  out.ind.sci = ( (ord.res.obs < sci.lower) | (ord.res.obs > sci.upper) )
  
  # construct data frames for quantities to be used
  #sci.p.in = data.frame(xin=ord.typ.res.sim[!out.ind.sci], yin=ord.res.obs[!out.ind.sci])
  sci.line = data.frame(lx=ord.typ.res.sim, ly.upp=sci.upper, ly.low=sci.lower)
  sci.p.out = data.frame(xout=ord.typ.res.sim[out.ind.sci], yout=ord.res.obs[out.ind.sci])
  
  
  ## ------------------------------------------------------------------------------------------
  ## Pointwise prediction interval
  ## ------------------------------------------------------------------------------------------
  # determine what quantiles to use for plotting:
  alpha = 1-envelope
  q.upp = envelope + alpha/2
  q.low = alpha/2
  quant = apply(ord.res.sim.mat[-(sim+1), ], 2, quantile, probs = c(q.low, q.upp))
  
  ## plotting options:
  min.x = min(ord.typ.res.sim)
  max.x = max(ord.typ.res.sim)
  min.y = min(quant[1, ], ord.res.obs)
  max.y = max(quant[2, ], ord.res.obs)
  
  #   # estiamte p.hat:
  out.ind.pw = ( (ord.res.obs < quant[1, ]) | (ord.res.obs > quant[2, ]) )
  # for plotting purposes, only those outside the pointwise CI but inside the simultaneous CI
  out.ind.pw2 = out.ind.pw==TRUE & out.ind.sci==FALSE
  #   nout = sum(out.ind)
  #   p.hat = nout/n.pts
  
  # construct data frames for quantities to be used
  epp.p.in = data.frame(xin=ord.typ.res.sim[!out.ind.pw], yin=ord.res.obs[!out.ind.pw])
  epp.line = data.frame(lx=ord.typ.res.sim, ly.low=quant[1, ], ly.upp=quant[2,])
  epp.p.out = data.frame(xout=ord.typ.res.sim[out.ind.pw2], yout=ord.res.obs[out.ind.pw2])
  
  
  # more information as legends
  info = data.frame(model=model, pv.D=pv.D, pv.P=pv.P)
  
  # return as a list
  return(list(epp.p.in = epp.p.in, epp.line = epp.line, epp.p.out = epp.p.out, 
              sci.line = sci.line, sci.p.out = sci.p.out, 
              info = info)
  )
}



#' @title Empirical Probability Plots of Negative Binomial Regression Models
#' 
#' @description This function makes empirical probability plots (EPP) based on the goodness-of-fit test results (an "gofv" object) from
#' \code{\link{nb.gof.v}}.
#' 
#' @param gofv.obj  an object of class "gofv" from the \code{\link{nb.gof.v}} output.
#' @param envelope  the level of prediction band when drawing EPP.
#' @param data.note  user-specified notes on lower right part of the plot.
#' 
#' @export
#' 
#' @usage EPPlot(gofv.obj, envelope=0.95, data.note=NULL)
#' 
#' @return A ggplot object of empirical probability plot.
#' 
#' @seealso \code{\link{nb.gof.v}} for simulated data examples.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
EPPlot = function(gofv.obj, envelope=0.95, data.note=NULL){
  
  eppdata = eppdata(gofv.obj, envelope=envelope)
  
  title = with(eppdata$info, paste("Test: ", model, ". Actual: ", data.note, sep=""))
  pv.D = ifelse(eppdata$info$pv.D >= 0.0001, 
                with(eppdata$info, paste("GOF p-value (Sq.Vert.Dist.) = ", as.numeric(formatC(pv.D, flag="#", digits=2)), sep="")),
                "GOF p-value (Sq.Vert.Dist.) < 0.0001")
  pv.P = ifelse(eppdata$info$pv.P >= 0.0001, 
                with(eppdata$info, paste("GOF p-value (Pear.Stat.) = ", as.numeric(formatC(pv.P, flag="#", digits=2)), sep="")),
                "GOF p-value (Pear.Stat.) < 0.0001")
  
  # begin plotting
  epp = ggplot(data=eppdata$epp.p.in, aes(x = xin, y = yin)) + 
    geom_point(data=eppdata$epp.p.in, aes(x = xin, y = yin), colour = "snow4", alpha = 1, size=1.5, shape=20) +
    geom_point(data=eppdata$epp.p.out, aes(x = xout, y = yout), colour = "Blue", alpha = 1, size=1.5, shape=2) + 
    geom_point(data=eppdata$sci.p.out, aes(x = xout, y = yout), colour = "Red", alpha = 1, size=1.5, shape=3) + 
    geom_line(data=eppdata$epp.line, aes(x = lx, y = ly.upp), linetype="dashed", colour="Blue") +
    geom_line(data=eppdata$epp.line, aes(x = lx, y = ly.low), linetype="dashed", colour="Blue") + 
    geom_line(data=eppdata$sci.line, aes(x = lx, y = ly.upp), linetype="solid", colour="Red") +
    geom_line(data=eppdata$sci.line, aes(x = lx, y = ly.low), linetype="solid", colour="Red") + 
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    theme_bw() +
    scale_x_continuous("") +
    scale_y_continuous("") +
    ggtitle(title) + 
    annotate("text", x = -Inf, y = Inf, hjust=0, vjust=2, label=pv.P, size = 3) +
    annotate("text", x = -Inf, y = Inf, hjust=0, vjust=4, label=pv.D, size = 3) +
    theme(plot.title = element_text(face="bold", size=11),
          axis.title.x = element_text(face="bold", size=8),
          axis.title.y = element_text(face="bold", size=8, angle=90),
          axis.text.x  = element_text(size=9),
          axis.text.y  = element_text(size=9),
          plot.margin=unit(x=c(0,0.2,0,0.2), units="cm")  # top, right, bottom, and left
    )
}
