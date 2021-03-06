
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
#' 
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
#' gof.nb2 = nb.gof.v(y, x, sim=999, model="NB2")
#' gof.nbp = nb.gof.v(y, x, sim=999, model="NBP", method="ML")
#' x2 = x[,2]
#' gof.poi = nb.gof.v(y, x2, sim=999, model="Poisson")
#' 
#' ## Empirical Probability Plots:
#' # pdf(file=file.path(path1, "eqk-nb-models.pdf"), width=8, height=4)
#' par(mfrow=c(1,2))
#' print(EPPlot(gof.nb2, envelope=0.95, data.note="Earthquake Dataset"))
#' print(EPPlot(gof.nbp, envelope=0.95, data.note="Earthquake Dataset"))
#' # dev.off()
#' 
NULL
