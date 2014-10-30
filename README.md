NBGOF
=====

**Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA sequencing Data**

----------------------------------------------------------------------------------------------------

This README file for the R package NBGOF is incomplete for now. We will update on a regular basis.

The source package can be downloaded [here](https://www.dropbox.com/s/lzpi5apn8el7may/NBGOF_0.1.6.tar.gz?dl=0)

### Install
Installation from this repository requires the `devtools` package pre-installed. Run the following R command to install `NBGOF`:

```S
devtools::install_github('NBGOF', 'gu-mi')
```

The NBGOF package implements goodness-of-fit (GOF) tests for negative binomial (NB) distributions and NB dispersion models, with applications in RNA-Seq data analysis. This package can be used to test the GOF of the NB2, NBP or Poisson **regression models**. It can also be used to test GOF for a variety of **NB dispersion models** in popular R/Bioconductor packages, including

* NBP dispersion model in the NBPSeq package (NBP)

* NBQ dispersion model in the NBPSeq package (NBQ)

* NB common dispersion model in the edgeR package (Common)

* NB genewise dispersion model in the edgeR package (Genewise)

* NB trended (non-parametric) dispersion model in the edgeR package (Trended)

* NB tagwise-common dispersion model in the edgeR package (Tagwise-Common)

* NB tagwise-trended dispersion model in the edgeR package (Tagwise-Trend)

To load the package into current R session, run

```S
library(NBGOF)
```
Two main functions for testing the adequacy of regression models and NB dispersion models are `nb.gof.v` and `nb.gof.m`, repectively. We provide two real datasets to illustrate the use of this package.

All results (figures/tables) in the submitted manuscript was produced using the following versions of R and its packages:

```{r}
R version 3.0.2 (2013-09-25)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      splines   parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_0.2.3    plyr_1.8.1      edgeR_3.2.4     limma_3.16.8    NBPSeq_0.2.9    qvalue_1.34.0  
 [7] NBGOF_0.1.4     ggplot2_0.9.3.1 doMC_1.3.3      iterators_1.0.6 foreach_1.4.1  

loaded via a namespace (and not attached):
 [1] codetools_0.2-8    colorspace_1.2-4   dichromat_2.0-0    digest_0.6.4       gtable_0.1.2      
 [6] labeling_0.2       MASS_7.3-29        munsell_0.4.2      numDeriv_2012.9-1  proto_0.3-10      
[11] RColorBrewer_1.0-5 Rcpp_0.11.0        reshape2_1.2.2     stringr_0.6.2      tcltk_3.0.2       
[16] tools_3.0.2       
```

**The current package version is 0.1.5, the beta version for next release including new dispersion models and minor bug fixes to version 0.1.4.**
