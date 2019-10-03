# pflamelet

This is an R package to build Persistence Flamelets, a functional summary of Multiscale Persistent Homology, and make inference on them. 

It implements all the methods introduced in [Persistence Flamelets: multiscale Persistent Homology for kernel density exploration](https://arxiv.org/abs/1709.07097) plus bonus features such as building bootstrap confidence bands on the Flamelet and cleaning it from noise. 

While it is not on [CRAN](https://cran.r-project.org) yet, you can install it as follows: 
``` r
# install.packages("devtools")
devtools::install_github("tulliapadellini/pflamelet")
```
