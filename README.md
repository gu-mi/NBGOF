NBGOF
=====

### Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA sequencing Data

`NBGOF` is an R package for implementing goodness-of-fit (GOF) tests for negative binomial (NB) distributions and NB dispersion models, with applications in RNA-Seq data analysis. This package can be used to test the GOF of the NB2, NBP or Poisson regression models. It can also be used to test GOF for a variety of NB dispersion models in popular R/Bioconductor packages, including

* NBP dispersion model in the NBPSeq package (NBP)

* NBQ dispersion model in the NBPSeq package (NBQ)

* NB common dispersion model in the edgeR package (Common)

* NB genewise dispersion model in the edgeR package (Genewise)

* NB trended (non-parametric) dispersion model in the edgeR package (Trended)

* NB tagwise-common dispersion model in the edgeR package (Tagwise-Common)

* NB tagwise-trended dispersion model in the edgeR package (Tagwise-Trend)

The methodologies are discussed in the manuscript **Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA Sequencing Data** (by Gu Mi, Yanming Di, and Daniel W. Schafer, PLOS ONE, 10(3)). The paper is freely available from this [link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119254). Functions in `NBGOF` have been used to generate all figures and tables displayed in the manuscript, plus some additional analysis tools for further investigations.

******

**The `NBGOF` source package can be downloaded [here](https://www.dropbox.com/s/4nnsx8ms9hhsjwz/NBGOF_0.2.1.tar.gz?dl=0).**

### Install
Installation from this repository requires the `devtools` package pre-installed. Run the following R command to install `NBGOF`:

```S
devtools::install_github('NBGOF', 'gu-mi')
```

To load the package into current R session, run

```S
library(NBGOF)
```

Two main functions for testing the adequacy of regression models and NB dispersion models are `nb.gof.v` and `nb.gof.m`, repectively. We provide two real datasets (`arab` and `earthquake`) to illustrate the use of this package.

******

We provide (Dropbox) links below to download R source codes and related supporting files for preparing the datasets and reproducing figures/tables in the manuscript. Some intermediate key results are also provided when necessary.

### Figures

* Figure 1 -- The mean-dispersion plot with six fitted dispersion models (common, NBP, NBQ, trended, tagwise-common and tagwise-trend) for the Arabidopsis RNA-Seq dataset (19,623 genes from three biological samples in the mock treatment group): [Download R files and results](https://www.dropbox.com/sh/x5quzc102xnjhqe/AADCi_1nBqaYwDlV7OIV0Bi_a?dl=0)
* Figure 2 -- Empirical probability plots with GOF test p-values for evaluating NB2 and NBP model fits on the earthquake dataset: [Download R files and results](https://www.dropbox.com/sh/kw4u8i2d3k4ie2h/AAAFmDwivEHSAuCpSrX2eKS-a?dl=0)
* Figure 3 -- Empirical probability plots and GOF p-values for testing NB2 (top row) and NBP (bottom) on four simulated datasets with sample size = 45: [Download R files and results](https://www.dropbox.com/sh/zen2m2cjsg1zh7v/AAAAXmj605u8hO1-dKcPFovva?dl=0)
* Figure 4 -- Uniform QQ plots of individual GOF test p-values for the Arabidopsis dataset (based on a random sample of 1,000 genes from six experimental units in two experimental groups): [Download R files and results](https://www.dropbox.com/sh/nvrbltddxnmkpjt/AABFYEO1SaEZD7LIm5Bf3UJGa?dl=0)
* SI Figure 1 -- Uniform QQ plots of individual GOF test p-values for the simulated ``NB2+noise'' dataset: [Download R files and results](https://www.dropbox.com/sh/lbipomac0s9t7gv/AAC7OhlfxyxL8KZelgOPzaVpa?dl=0)


### Tables
* Table 1 -- Type I error rates for 0.05-level NB2 and NBP GOF tests: [Download R files](https://www.dropbox.com/sh/131q2kw0skoy3du/AABYvOBGEPodEjOqB9cUbNlia?dl=0)
* Table 2 -- Rejection rates for 0.05-level NB2 and NBP GOF tests: [Download R files](https://www.dropbox.com/sh/iidjc6x00ufzl0m/AADzJzZwneq203TFTFn9ipTUa?dl=0)

******

All results (figures/tables) in the submitted manuscript was produced using the following versions of R and its packages:

```{r}
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11 (El Capitan)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] NBGOF_0.2.0     ggplot2_1.0.1   doMC_1.3.3      iterators_1.0.7 foreach_1.4.2  

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.1       NBPSeq_0.3.0      magrittr_1.5      edgeR_3.10.5      splines_3.2.2     MASS_7.3-44      
 [7] munsell_0.4.2     colorspace_1.2-6  stringr_1.0.0     plyr_1.8.3        tools_3.2.2       gtable_0.1.2     
[13] digest_0.6.8      numDeriv_2014.2-1 reshape2_1.4.1    codetools_0.2-14  qvalue_2.0.0      labeling_0.3     
[19] limma_3.24.15     stringi_0.5-5     compiler_3.2.2    scales_0.3.0      proto_0.3-10     
```

******

If you have any questions, please do not hesitate to email the repository maintainer (Gu Mi) at neo.migu@gmail.com. Thank you for your interests in our research work.
