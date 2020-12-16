#' pflamelet: Persistence Flamelets
#'
#'
#' @docType package
#' @name pflamelet
#' @description Computes the Persistence Flamelets, a statistical tool for exploring the
#' Topological Invariants of Scale-Space families introduced in Padellini and Brutti (2017)
#' "Persistence Flamelets: multiscale Persistent Homology for kernel density exploration" <arXiv:1709.07097>.Flamelets can be built from either sub/super level sets of arbitrary functions or
#' from a precomputed list of Persistence Diagrams.
#' In addition, this package provides functions to compute confidence bands for a Flamelet via bootstrap, assess significance of each feature,
#' clean a Flamelet from topological noise, perform a two sample permutation test for groups of Flamelets and it also implements a topological heuristic for bandwidth selection for kernel density estimators.
#' @importFrom stats quantile
#' @importFrom TDA gridDiag kde landscape silhouette distFct bootstrapDiagram
#' @importFrom pbapply pblapply pbreplicate
#' @importFrom plotly plot_ly add_surface layout add_trace  %>%
#' @importFrom abind abind
#' @importFrom graphics image
#' @importFrom viridis viridis
NULL
