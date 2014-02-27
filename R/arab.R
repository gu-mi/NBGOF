
#' @title The Arabidopsis RNA-Seq Dataset
#'
#' @description An RNA-Seq dataset from a pilot study of the defense response of
#' Arabidopsis to infection by bacteria. We performed RNA-Seq
#' experiments on three independent biological samples from each of
#' the two treatment groups.  The matrix contains the frequencies of
#' RNA-Seq reads mapped to genes in a reference database. Rows
#' correspond to genes and columns correspond to independent
#' biological samples.
#' 
#' @details   We challenged leaves of Arabidopsis with the defense-eliciting
#' \emph{\eqn{\Delta}hrcC} mutant of \emph{Pseudomonas syringae} pathovar
#' \emph{tomato} DC3000.  We also infiltrated leaves of Arabidopsis with 10mM
#' MgCl2 as a mock inoculation. RNA was isolated 7 hours after inoculation, 
#' enriched for mRNA and prepared for RNA-Seq. We sequenced one replicate per 
#' channel on the Illumina Genome Analyzer (http://www.illumina.com). 
#' The length of the RNA-Seq reads can vary in length depending on user preference 
#' and the sequencing instrument. The dataset used here are derived from a 36-cycle 
#' sequencing reaction, that we trimmed to 25mers.  We used an in-house computational 
#' pipeline to process, align, and assign RNA-Seq reads to genes according to a 
#' reference database we developed for Arabidopsis.
#' 
#' @usage data(arab)
#' @format A 26222 by 6 matrix of RNA-Seq read frequencies.
#' @name arab
#'
#' @docType data
#'
#' @references 
#' Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
#' Model for Assessing Differential Gene Expression from RNA-Seq", \emph{Statistical 
#' Applications in Genetics and Molecular Biology}, 10 (1).
#'
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
#' @keywords datasets
#'
#' @examples
#' 
#' ## Load the dataset into R session:
#' library(NBGOF)
#' data(arab)
#' 
#' ## simulation set-up:
#' seed = 539
#' sim = 999
#' conf.env = 0.95
#' nc = detectCores() - 1
#' set.seed(seed)
#' 
#' ## subset read counts to exclude all zero rows within each treatment:
#' good.arab = arab[(rowSums(arab[ ,1:3]) != 0 & rowSums(arab[ ,4:6]) != 0), ]
#' m = 500    # number of genes retained
#' # consider a single group with 3 replicates
#' samp.idx = sample(1:dim(good.arab)[1], m)
#' arab.sub = good.arab[samp.idx,1:3] 
#' lib.sizes = colSums(arab.sub)
#' y = arab.sub
#' 
#' ## model matrix for arab.sub:
#' x = as.matrix(rep(1,3))
#' 
#' ## GOF tests for different dispersion models:
#' fnbp.arab = nb.gof.m(counts=y, x=x, sim=sim, model="NBP", seed=1, ncores=nc)
#' fnbq.arab = nb.gof.m(counts=y, x=x, sim=sim, model="NBQ", seed=1, ncores=nc)
#' fcom.arab = nb.gof.m(counts=y, x=x, sim=sim, model="Common", seed=1, ncores=nc)
#' fgen.arab = nb.gof.m(counts=y, x=x, sim=sim, model="Genewise", seed=1, ncores=nc)
#' ftgc.arab = nb.gof.m(counts=y, x=x, sim=sim, model="Tagwise-Common", prior.df=10, seed=1, ncores=nc)
#' ftgt.arab = nb.gof.m(counts=y, x=x, sim=sim, model="Tagwise-Trend", prior.df=10, min.n=100, seed=1, ncores=nc)
#' ftrd.arab = nb.gof.m(counts=y, x=x, sim=sim, model="Trended", min.n=100, seed=1, ncores=nc)
#' 
#' ## summarize the GOF test results:
#' summary(fnbp.arab, conf.env=conf.env, data.note="arab")
#' summary(fnbq.arab, conf.env=conf.env, data.note="arab")
#' summary(fcom.arab, conf.env=conf.env, data.note="arab")
#' summary(fgen.arab, conf.env=conf.env, data.note="arab")
#' summary(ftgc.arab, conf.env=conf.env, data.note="arab")
#' summary(ftgt.arab, conf.env=conf.env, data.note="arab")
#' summary(ftrd.arab, conf.env=conf.env, data.note="arab")
#' 
#' ## mean-dispersion plot with some of the fitted dispersion models:
#' # model.vec=c("Common", "NBP", "NBQ", "Trended", "Tagwise-Common", "Tagwise-Trend")
#' # pl.arab = MDPlot(model.vec, y, x, title="Mean-Dispersion Plot with Fitted Dispersion Models (Arabidopsis Data)")
#' # pdf(file="arab-md-plot-subset.pdf", width=7, height=7)
#' # print(pl.arab)
#' # dev.off()
#' 
#' ## diagnostic plot (histogram or uniform QQ plot of GOF p-values):
#' # dp1 = diagPlot(fnbq.arab, type="hist", binwidth=0.1)
#' # print(dp1)
#  # dp2 = diagPlot(fnbq.arab, type="qq", envelope=0.95)
#  # print(dp2)
#' 
NULL

