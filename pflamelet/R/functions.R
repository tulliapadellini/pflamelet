# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Persistence Flamelet Function
#'
#' Computes the Persistence Flamelets from a list of Persistence Diagrams. If X is a matrix or a dataframe the function computes the flamelets of the
#' kernel density estimator computed on X.
#'
#'@param X   a list of persistence diagrams representing the scale-space family at different resolution.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the corresponding KDE with the bandwidth as scale parameter.
#'@param base.type a string specifying whether the Flamelet is built from Persistence Landscapes ("landscape") or Persistence Silhouettes ("silhouettes")
#'@param base.param the order k of the Flamelet (if base.type=="landscape") or the power p of the Flamelet (if base.type=="silhouette")
#'@param dimension the topological dimension of the flamelet (0 for connected components, 1 for loops, ...)
#'@param tseq a vector of values at which the Flamelet function is evaluated for a fixed scale level
#'@param h.grid vector of bandwidths for the KDE, representing the scale parameter of the Flamelet
#'@param lim 2-by-d matrix, where the i-th column contains the range of the grid over which KDE is computed for the i-th variable.
#'@param by a scalar (or a vector if different values are selected for each dimension)
#'             specifying spaces between elements on the grid
#'@param scale a logical indicating whether or not the Persistence Diagrams have to be scaled to be in the same range (needed only for visualization purposes)
#'@examples
#'
#'## library(TDA)
#'## xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'## Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'## lim = cbind(Xlim, Ylim)
#'## foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#'## base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
# ##                              tseq = seq(0, .75, length.out = 500))
#'@export
build.flamelet = function(X, base.type = "landscape", base.param = 1, dimension=1,
                          tseq, diag.fun = distFct, h.grid =NULL, lim = NULL, by=NULL,
                          scale = TRUE, precomputed.diagram = FALSE, sublevel = TRUE, info.message = FALSE){

  if(!is.list(X)){
    f = function(h) {
      diagr = gridDiag(X, kde, lim = lim, by = by, sublevel = FALSE,
                       location = FALSE, printProgress = FALSE, h = h, maxdimension = dimension)$diagram
      if(scale) diagr[,2:3] = diagr[,2:3]/diagr[1,3]
      return(diagr)
    }

    if(info.message) cat("Computing the Persistence Diagrams - hold on, this is the longest step...")
    diagList = lapply(h.grid, f)
    if(info.message) cat(" and now the Persistence Landscapes")

  } else {
    if(precomputed.diagram) diagListh = X
    else{

      f = function(x) {
        diagr = gridDiag(x, diag.fun, lim = lim, by = by, sublevel = sublevel,
                         location = FALSE, printProgress = FALSE, maxdimension = dimension)$diagram
        return(diagr)
      }

      if(info.message) cat("Computing the Persistence Diagrams - hold on, this is the longest step...")
      diagList = lapply(X, f)
      if(info.message) cat(" and now the Persistence Landscapes")

    }
  }

  if(base.type == "landscape") flamelet = sapply(diagList, landscape, dimension = dimension, KK = base.param, tseq = tseq)
  if(base.type == "silhouette") flamelet = sapply(diagList, silhouette, dimension = dimension, p = base.param, tseq = tseq)
  return(flamelet)

}



#' Bootstrap Band for Persistence Flamelet
#'
#' Computes a bootstrap band for the Persistence Flamelet from a list of Persistence Diagrams. If X is a matrix or a dataframe the function computes the flamelets of the
#' kernel density estimator computed on X.
#'
#'@param X   a list of persistence diagrams representing the scale-space family at different resolution.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the corresponding KDE with the bandwidth as scale parameter.
#'@param base.type a string specifying whether the Flamelet is built from Persistence Landscapes ("landscape") or Persistence Silhouettes ("silhouettes")
#'@param base.param the order k of the Flamelet (if base.type=="landscape") or the power p of the Flamelet (if base.type=="silhouette")
#'@param dimension the topological dimension of the flamelet (0 for connected components, 1 for loops, ...)
#'@param tseq a vector of values at which the Flamelet function is evaluated for a fixed scale level
#'@param h.grid vector of bandwidths for the KDE, representing the scale parameter of the Flamelet
#'@param lim 2-by-d matrix, where the i-th column contains the range of the grid over which KDE is computed for the i-th variable.
#'@param by a scalar (or a vector if different values are selected for each dimension)
#'             specifying spaces between elements on the grid
#'@param scale a logical indicating whether or not the Persistence Diagrams have to be scaled to be in the same range (needed only for visualization purposes)
#'@examples
#'
#'## library(TDA)
#'## xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'## Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'## lim = cbind(Xlim, Ylim)
#'## foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#'## base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
# ##                              tseq = seq(0, .75, length.out = 500))
#'@export
flamelet.band = function(X, B, alpha, base.type = "landscape", base.param = 1, dimension=1,
                         tseq, diag.fun = distFct, h.grid =NULL, lim = NULL, by=NULL,
                         sublevel = TRUE){
  # X: list of data
  # diag.build: function to build the diagrams
  # B: number of bootstrap rep

  # if X is a matrix


  if(is.list(X)){

    Xbase = build.flamelet(X, base.type = base.type, base.param = base.param, dimension=,
                        tseq=tseq, diag.fun = diag.fun, h.grid = h.grid, lim = lim, by= by,
                        scale = FALSE, precomputed.diagram = FALSE, sublevel = sublevel, info.message = FALSE)
  f = function(){
    subset.idx = sample(1:nrow(X[[1]]), size = nrow(X[[1]]), replace = T)
    Xs = lapply(X, function(x) x[subset.idx, ])
    Xf = build.flamelet(Xs, base.type = base.type, base.param = base.param, dimension=,
                                tseq=tseq, diag.fun = diag.fun, h.grid = h.grid, lim = lim, by= by,
                                scale = FALSE, precomputed.diagram = FALSE, sublevel = sublevel, info.message = FALSE)
    return(max(abs(Xf-Xbase)))

  }
  } else{
    Xbase = build.flamelet(X, base.type = base.type, base.param = base.param, dimension=,
                        tseq=tseq, diag.fun = diag.fun, h.grid = h.grid, lim = lim, by= by,
                        scale = FALSE, precomputed.diagram = FALSE, sublevel = sublevel, info.message = FALSE)

    f = function(){
      subset.idx = sample(1:nrow(X), size = nrow(X), replace = T)
      Xs = X[subset.idx, ]
      Xf = build.flamelet(Xs, base.type = base.type, base.param = base.param, dimension=,
                          tseq=tseq, diag.fun = diag.fun, h.grid = h.grid, lim = lim, by= by,
                          scale = FALSE, precomputed.diagram = FALSE, sublevel = sublevel, info.message = FALSE)
      return(max(abs(Xf-Xbase)))

    }



  }

  bvect = pbapply::pbreplicate(B, f())

  return(quantile(bvect, alpha))


}


# if (plot) {
#   three.d = plot_ly(x = h.grid, y = tseq,  z = landSurf, colors = viridis(100, option = "B", end = .87)) %>%
#     add_surface() %>%
#     layout(
#       scene = list(
#         xaxis = list(title = "Bandwidth"),
#         yaxis = list(title = "(Birth + Death)/2"),
#         zaxis = list(title = "Persistence"))
#     )
