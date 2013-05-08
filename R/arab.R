#' The Arabidopsis RNA-Seq Data Set
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
#' @references Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The NBP Negative Binomial
#' Model for Assessing Differential Gene Expression from RNA-Seq", Statistical Applications in
#' Genetics and Molecular Biology, 10 (1).
#'
#' @keywords datasets
#'
#' @examples
#' ## Load the dataset into R session:
#' data(arab)
#' 
#' ## simulation set-up
#' seed = 315826
#' sim = 999
#' conf.env = 0.95
#' set.seed(seed)
#' 
#' ## subset read counts to exclude all zero rows within each treatment
#' good.arab = arab[(rowSums(arab[ ,1:3]) != 0 & rowSums(arab[ ,4:6]) != 0), ]
#' m = 500    # number of genes subsetted
#' # consider a single group with 3 replicates
#' samp.idx = sample(1:dim(good.arab)[1], m)  # sample m genes
#' arab.sub = good.arab[samp.idx,1:3] 
#' lib.sizes = colSums(arab.sub)
#' 
#' # model matrix for arab.sub
#' x = as.matrix(rep(1,3))
#' x
#' 
#' # GOF tests for different dispersion models
#' fnbp.nbp = nb_gof_m(counts=y, x=x, sim=sim, model.fit="NBP")
#' fcom.nbp = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-common")
#' fgen.nbp = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-genewise")
#' ftag.nbp = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-tagwise")
#' ftrd.nbp = nb_gof_m(counts=y, x=x, sim=sim, model.fit="edgeR-trended", min.n=100)
#' 
#' # summarize the GOF test results
#' summary(fnbp.nbp, conf.env=conf.env, data.note="arab")
#' summary(fcom.nbp, conf.env=conf.env, data.note="arab")
#' summary(fgen.nbp, conf.env=conf.env, data.note="arab")
#' summary(ftag.nbp, conf.env=conf.env, data.note="arab")
#' summary(ftrd.nbp, conf.env=conf.env, data.note="arab")
#' 
NULL