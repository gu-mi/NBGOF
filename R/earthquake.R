
#' @title The Earthquake Dataset
#' 
#' @description This dataset contains the frequencies of all earthquakes of a given magnitude
#' (reported to one decimal place) for magnitudes from 4.5 to 9.1, that occurred between 
#' January 1, 1964 to December 31, 2012.
#' 
#' @details We use this dataset as a real data example to illustrate goodness-of-fit tests
#' of NB and Poisson regression models (univariate response).
#' 
#' @usage data(earthquake)
#' @format A 45 by 2 matrix, with column names "Magnitude" and "Frequency"
#' @name earthquake
#'
#' @docType data
#'
#' @references Composite Earthquake Catalog, Advanced National Seismic System, Northern 
#' California Earthquake Data Center (NCEDC), \url{http://quake.geo.berkeley.edu/cnss/}.
#' 
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#'
#' @keywords datasets
#'
#' @examples
#' ## Load the dataset into R session:
#' library(NBGOF)
#' data(earthquake)
#' 
#' ## basic descriptions of the dataset:
#' head(earthquake)
#' range(earthquake$Magnitude)  # 4.5  9.1
#' range(earthquake$Frequency)  # 1  33280
#' 
#' ## GOF test of NB2, NBP and Poisson models:
#' y = earthquake$Frequency
#' x = as.matrix(cbind(rep(1,length(y)), earthquake$Magnitude))
#' gf.nb2 = nb_gof_v(y, x, sim=999, model="NB2")
#' gf.nbp = nb_gof_v(y, x, sim=999, model="NBP", est.method="ML")
#' x2 = x[,2]
#' gf.poi = nb_gof_v(y, x2, sim=999, model="Poisson")
#' 
#' ## empirical probability plots:
#' # pdf(file=file.path(path1, "eqk-nb-models-12.pdf"), width=8, height=4)
#' par(mfrow=c(1,3))
#' plot(gf.nb2, pch=".", cex=5,  leg.cex = 0.9, data.note="Earthquake Dataset")
#' plot(gf.nbp, pch=".", cex=5,  leg.cex = 0.9, data.note="Earthquake Dataset")
#' plot(gf.poi, pch=".", cex=5,  leg.cex = 0.9, data.note="Earthquake Dataset")
#' # dev.off()
NULL

