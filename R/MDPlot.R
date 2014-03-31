
#' @title Prepare Quantities for Mean-Dispersion Plots Using ggplot2
#' 
#' @description This function prepares quantities for making mean-dispersion plots with curves of fitted negative binomial dispersion models. The result (a data frame) will be passed to the \code{\link{MDPlot}} function.
#' 
#' @param counts an m-by-n count matrix of non-negative integers. For a typical RNA-Seq experiment, this is the read counts with 
#' m genes and n samples
#' @param x an n-by-p design matrix
#' @param model a string of characters specifying the negative binomial dispersion model used to fit the data. Currently supported
#' dispersion models include "NBP", "NBQ", "NBS", "STEP", "Common", "Tagwise-Common", "Tagwise-Trend" and "Trended"
#' 
#' @return a data frame for each specified model.
#' 
#' @keywords internal
#' 
mddata = function(counts, x, model = NULL){
  
  m = dim(counts)[1]
  n = dim(counts)[2]
  
  # naive estimations of relative frequency and phi for the basic scatter plot
  # it makes sense to have only one group (design matrix x is a column: intercept-only model)
  
  mu.mom = edgeR:::expandAsMatrix(rowMeans(counts), dim=c(m,n))
  re.freq = (mu.mom / (matrix(1, m, 1) %*% matrix(colSums(counts), 1, n)))[ ,1] 
  #qt.re.freq = quantile(re.freq, c(0.001, 0.999))
  phi.mom = (rowSums((counts - mu.mom)^2) - rowSums(mu.mom))/rowSums(mu.mom^2)  # may have NaN
  #id = (phi.mom > 0 & !is.nan(phi.mom) & re.freq > qt.re.freq[1] & re.freq < qt.re.freq[2])  
  id = (phi.mom > 0 & !is.nan(phi.mom))
  # may discard some phi.hat here not in the plotting: this "id" is used "globally" to subset
  n.pt = sum(id)
  
  #### -----------------------------------------------------------------
  if (is.null(model)){
    mu.mom = edgeR:::expandAsMatrix(rowMeans(counts), dim=c(m,n))
    re.freq = (mu.mom / (matrix(1, m, 1) %*% matrix(colSums(counts), 1, n)))[ ,1] 
    #qt.re.freq = quantile(re.freq, 0.001, 0.999)
    phi.mom = (rowSums((counts - mu.mom)^2) - rowSums(mu.mom))/rowSums(mu.mom^2)  # may have NaN
    #id = (phi.mom > 0 & !is.nan(phi.mom) & re.freq > qt.re.freq[1] & re.freq < qt.re.freq[2])   
    id = (phi.mom > 0 & !is.nan(phi.mom))
    # may discard some phi.hat here not in the plotting
    coords.scatter = data.frame(x0=re.freq[id], y0=phi.mom[id])   # for scatter plot
    return(coords.scatter)
  }
  
  #### -----------------------------------------------------------------
  if (model == "NBP"){
    nb.data = prepare.nb.data(counts = counts, lib.sizes = colSums(counts),norm.factors = rep(1, n)) 
    nbp.disp.z = disp.nbp(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbp.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NBP", method = "MAPL")    
    phi.nbp = nbp.disp$estimates
    pi.hat = nbp.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbp[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("NBP", m), coords)))
  }
  
  #### -----------------------------------------------------------------
  if (model == "NBQ"){
    nb.data = prepare.nb.data(counts = counts, lib.sizes = colSums(counts), norm.factors = rep(1, n))  
    nbq.disp.z = disp.nbq(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbq.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NBQ", method = "MAPL")
    phi.nbq = nbq.disp$estimates
    pi.hat = nbq.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbq[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("NBQ", m), coords)))  
  }
  
  #### -----------------------------------------------------------------
  if (model == "NBS"){
    nb.data = prepare.nb.data(counts = counts, lib.sizes = colSums(counts), norm.factors = rep(1, n))  
    nbs.disp.z = disp.nbs(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbs.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NBS", method = "MAPL")
    phi.nbs = nbs.disp$estimates
    pi.hat = nbs.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbs[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("NBS", m), coords)))  
  }
  
  #### -----------------------------------------------------------------
  if (model == "STEP"){
    nb.data = prepare.nb.data(counts = counts, lib.sizes = colSums(counts), norm.factors = rep(1, n))  
    nbstep.disp.z = disp.step(counts=counts, eff.lib.sizes=nb.data$eff.lib.sizes, x=x)   # to get "z"
    nbstep.disp = estimate.dispersion(nb.data = nb.data, x = x, model = "NB2", method = "MAPL")
    phi.nbstep = nbstep.disp$estimates
    pi.hat = nbstep.disp.z$pi.pre[,1]
    coords = data.frame(x=pi.hat, y=phi.nbstep[ ,1])
    return(as.data.frame(cbind(Dispersion.Model = rep("STEP", m), coords)))  
  }
  
  #### -----------------------------------------------------------------
  if (model == "Common"){
    y.dge = DGEList(counts)
    e.com = estimateGLMCommonDisp(y.dge, x)
    phi = rep(e.com$common.dispersion, m)
    rel.freq = exp(e.com$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)   
    return(as.data.frame(cbind(Dispersion.Model = rep("Common", m), coords)))
  }
  
  #### -----------------------------------------------------------------
  if (model == "Tagwise-Common"){
    
    y.dge = DGEList(counts)
    e.com = estimateGLMCommonDisp(y.dge, x)
    e.tgc = estimateGLMTagwiseDisp(y.dge, design=x, dispersion=e.com$common.dispersion, trend=FALSE)
    phi = e.tgc$tagwise.dispersion
    rel.freq = exp(e.tgc$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)   
    return(as.data.frame(cbind(Dispersion.Model = rep("Tagwise-Common", m), coords)))
  }
  
  #### -----------------------------------------------------------------      
  if (model == "Tagwise-Trend"){
    
    y.dge = DGEList(counts)
    e.trd = estimateGLMTrendedDisp(y.dge, x)
    e.tgt = estimateGLMTagwiseDisp(y.dge, design=x, dispersion=e.trd$trended.dispersion, trend=TRUE)
    phi = e.tgt$tagwise.dispersion
    rel.freq = exp(e.tgt$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)    
    return(as.data.frame(cbind(Dispersion.Model = rep("Tagwise-Trend", m), coords))) 
  }
  
  #### -----------------------------------------------------------------
  if (model == "Trended"){ 
    trd = dispBinTrend(counts)
    # names(trd)  # "AveLogCPM"      "dispersion"     "bin.AveLogCPM"  "bin.dispersion"
    phi = trd$dispersion  # phi without noise
    rel.freq = exp(trd$AveLogCPM*log(2) - log(1e+06))
    # for plotting purpose:
    id2 = order(rel.freq)
    rel.freq.ord = rel.freq[id2]
    phi.ord = phi[id2]
    coords = data.frame(x=rel.freq.ord, y=phi.ord)
    return(as.data.frame(cbind(Dispersion.Model = rep("Trended", m), coords)))
  }
}


#' @title Mean-Dispersion Plots with Curves of Fitted Negative Binomial Dispersion Models
#' 
#' @description This function makes mean-dispersion plots with curves of fitted negative binomial dispersion models. We currently consider
#' only one group and fit an intercept-only model.
#' 
#' @param model.vec a character vector specifying the names of the NB dispersion models. Currently supported include 
#' "NBP", "NBQ", "Common", "Tagwise-Common", "Tagwise-Trend" and "Trended"
#' @param counts an m-by-n count matrix of non-negative integers. For a typical RNA-Seq experiment, this is the read counts with 
#' m genes and n samples
#' @param x an n-by-p design matrix
#' @param title title of the plot
#' 
#' @export
#' 
#' @usage MDPlot(model.vec, counts, x, title=NULL)
#' 
#' @return A ggplot object of mean-dispersion plot with fitted NB dispersion model curves.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
#' @examples
#' 
#' # load package
#' library(NBGOF)
#' 
#' # load data and subset for illustration
#' data(arab)
#' counts = arab[1:5000,1:3]
#' x = matrix(c(1,1,1), nr=3)
#' 
#' # specify what models to plot
#' model.vec=c("Common", "NBP", "NBQ", "Trended", "Tagwise-Common", "Tagwise-Trend")
#' 
#' # begin plotting
#' pl.arab = MDPlot(model.vec, counts, x, title="Mean-Dispersion Plot with Fitted Dispersion Models (Arabidopsis Data)")
#' print(pl.arab)
#' 
MDPlot = function(model.vec, counts, x, title=NULL, data.note=NULL){
  
  k = length(model.vec)
  result.lst = rep( list(NA), k) 
  for (i in seq_len(k)){
    result.lst[[i]] = mddata(counts, x, model = model.vec[i])
  }
  df.long = do.call(rbind.data.frame, result.lst)
  
  md = ggplot(data = mddata(counts, x, model=NULL), aes(x = x0, y = y0)) + 
    geom_point(data = mddata(counts, x, model=NULL), aes(x = x0, y = y0), alpha=I(.6), size=I(1)) + 
    scale_x_log10("Estimated Relative Frequency", breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_y_log10("Estimated NB Dispersion Parameter", breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    geom_line(data = df.long, aes(x=x, y=y, colour=factor(Dispersion.Model), size=factor(Dispersion.Model),
                                  linetype=factor(Dispersion.Model), alpha=factor(Dispersion.Model))) +
    scale_size_manual(values = c(1.5, 1.5, 1.5, 1.5, 0.5, 0.5)) + 
    scale_linetype_manual(values = seq(1, length(unique(df.long$Dispersion.Model)))) + 
    scale_alpha_manual(values = c(1, 1, 1, 1, 0.4, 0.4)) + 
    theme_bw() +
    ggtitle(title) + 
    annotate("text", x = 8e-7, y = 1e-5, label=data.note, size = 5) +
    theme(plot.title = element_text(face="bold", size=16),
          axis.title.x = element_text(face="bold", size=16),
          axis.title.y = element_text(face="bold", size=16, angle=90),
          axis.text.x  = element_text(size=12),
          axis.text.y  = element_text(size=12),
          legend.key.width = unit(3, "line"),
          legend.justification=c(1,0), 
          legend.position=c(1,0),
          legend.title = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size=14)
    )
}
