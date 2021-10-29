# SSLB R Package

This repository contains the R package SSLB, which implements the method Spike-and-Slab Lasso Biclustering (Moran et al, 2020).

## Installation

This package requires the R packages Rcpp (>= 1.0.5), RcppArmadillo and MASS.

To install the latest version of the package from Github, copy and paste the following code in your R GUI.
```
library(devtools)
install_github("gemoran/SSLB")
```

To install the package from the source code `SSLB.zip` in the supplementary material to Moran et al. (2021):

1. Expand the zip file `SSLB.zip`. The folder `SSLB` will result. 
2. In your R GUI, set your working directory to the folder `SSLB`.
3. Run the following R code:
```
library(devtools)
build()
install()
```

## Reference

Moran, G. E., Rockova, V., & George, E. I. (2021). Spike-and-slab lasso biclustering. The Annals of Applied Statistics, 15(1), 148-173.

DOI: 10.1214/20-AOAS1385