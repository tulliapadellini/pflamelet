# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Persistence Flamelet Function
#'
#' Computes the Persistence Flamelets from a list of Persistence Diagrams. If X is a matrix or a dataframe the function computes the Flamelets of the sub/superlevel set of
#' an arbitrary function computed on X.
#'
#'@param X   a list of persistence diagrams representing the scale-space family at different resolution.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the \code{diag.fun}.
#'@param base.type a string specifying whether the Flamelet is built from Persistence Landscapes ("landscape") or Persistence Silhouettes ("silhouettes")
#'@param base.param the order k of the Flamelet (if base.type=="landscape") or the power p of the Flamelet (if base.type=="silhouette")
#'@param dimension the topological dimension of the flamelet (0 for connected components, 1 for loops, ...)
#'@param tseq a vector of values at which the Flamelet function is evaluated for a fixed scale level
#'@param diag.fun the function whose sub/super-level set define the persistent homology groups. Corresponds to the argument \code{FUN} of the \code{gridDiag} function in the TDA package.
#'@param sublevel a logical indicating whether the Persistent Homology should be computed on sub or superlevel set of the function given as \code{diag.fun}
#'@param h.grid vector of bandwidths for the KDE, representing the scale parameter of the Flamelet
#'@param lim 2-by-d matrix, where the i-th column contains the range of the grid over which the function specified in \code{diag.fun} is computed for the i-th variable.
#'@param by a scalar (or a vector if different values are selected for each dimension)
#'             specifying spaces between elements on the grid whose outernmost element are defined by \code{lim}.
#'@param scale a logical indicating whether or not the Persistence Diagrams have to be scaled to be in the same range (needed only for visualization purposes)
#'@param precomputed.diagram a logical indicating whether the Persistence Diagrams have to be computed or are given as input in the form of a list of diagrams
#'@param info.message a logical denoting if progress messages should be printed
#'@param band a logical indicating whether or not each Persistence Diagram should be cleaned by means of Bootstrap Bands (only possible when the Persistence Diagrams are computed within the function and are not given as an input)
#'@param B number of bootstrap repetitions needed to compute the confidence band over Persistence Diagrams
#'@param alpha the confidence level of the bootstrap confidence bands
#'
#'
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
                           tseq, diag.fun = distFct, sublevel = TRUE, h.grid =NULL, lim = NULL, by=NULL,
                           scale = TRUE, precomputed.diagram = FALSE, info.message = FALSE, band = FALSE, B = 10, alpha = 0.95){

  if(!is.list(X)){
    f = function(h) {
      diagr = gridDiag(X, kde, lim = lim, by = by, sublevel = FALSE,
                       location = FALSE, printProgress = FALSE, h = h, maxdimension = dimension)$diagram

      if(nrow(diagr)>0){
        if(band) {
          cc <- bootstrapDiagram(X, kde, lim = lim, by = by, sublevel = FALSE, B = B,
                                 alpha = 1-alpha/length(h.grid), dimension = dimension, printProgress = FALSE, h = h)


          idx = apply(as.matrix(diagr[,2:3]), 1, diff)>2*cc
          diagr = diagr[idx,]
          diagr = matrix(diagr, ncol = 3)

        }
      }

      if(scale && nrow(diagr)>0) diagr[,2:3] = diagr[,2:3]/diagr[1,3]

      return(diagr)
    }

    if(info.message) cat("Computing the Persistence Diagrams - hold on, this is the longest step...")
    diagList = pbapply::pblapply(h.grid, f)
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
      diagList = pbapply::pblapply(X, f)
      if(info.message) cat(" and now the Persistence Landscapes")

    }
  }

  if(base.type == "landscape") flamelet = sapply(diagList, landscape, dimension = dimension, KK = base.param, tseq = tseq)
  if(base.type == "silhouette") flamelet = sapply(diagList, silhouette, dimension = dimension, p = base.param, tseq = tseq)
  return(flamelet)

}


#' Bootstrap Band for Persistence Flamelet
#'
#' Computes a bootstrap band for the Persistence Flamelet from a list of Persistence Diagrams. If X is a matrix or a dataframe the function computes the Flamelets of the sub/superlevel set of
#' an arbitrary function computed on X.
#'
#'@param X   a list of persistence diagrams representing the scale-space family at different resolution.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the corresponding \code{diag.fun}.
#'@param B number of bootstrap repetitions needed to compute the confidence band over Persistence Diagrams
#'@param alpha the confidence level of the bootstrap confidence bands
#'@param base.type a string specifying whether the Flamelet is built from Persistence Landscapes ("landscape") or Persistence Silhouettes ("silhouettes")
#'@param base.param the order k of the Flamelet (if base.type=="landscape") or the power p of the Flamelet (if base.type=="silhouette")
#'@param dimension the topological dimension of the flamelet (0 for connected components, 1 for loops, ...)
#'@param tseq a vector of values at which the Flamelet function is evaluated for a fixed scale level
#'@param diag.fun the function whose sub/super-level set define the persistent homology groups. Corresponds to the argument \code{FUN} of the \code{gridDiag} function in the TDA package.
#'@param sublevel a logical indicating whether the Persistent Homology should be computed on sub or superlevel set of the function given as \code{diag.fun}
#'@param h.grid vector of bandwidths for the KDE, representing the scale parameter of the Flamelet
#'@param lim 2-by-d matrix, where the i-th column contains the range of the grid over which the function specified in \code{diag.fun} is computed for the i-th variable.
#'@param by a scalar (or a vector if different values are selected for each dimension)
#'             specifying spaces between elements on the grid whose outernmost element are defined by \code{lim}.
#'#'@examples
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
                         tseq, diag.fun = distFct, sublevel = TRUE, h.grid =NULL, lim = NULL, by=NULL){


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



#' Plot Persistence Flamelet
#'
#' Plots
#'
#'@param X   a list of persistence diagrams representing the scale-space family at different resolution.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the \code{diag.fun}.
#'@param tseq a vector of values at which the Flamelet function is evaluated for a fixed scale level
#'@param scale.param a vector corresponding to the values of the scale parameter at which the Flamelet \code{X} has been evaluated
#'
#'
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
flamelet.plot = function(X, scale.param, tseq, flat,  scale.name ="Bandwidth"){


  if(flat){
  image(y = scale.param , x = tseq, z= X,
        col = viridis(100, option = "B"), ylab = scale.name, xlab = "(Birth + Death)/2")

  }else{
    out.plot = plot_ly(x = scale.param, y = tseq,  z = X, colors = viridis(100, option = "B", end = .87)) %>%
      add_surface() %>%
      layout(
        scene = list(
          xaxis = list(title = scale.name),
          yaxis = list(title = "(Birth + Death)/2"),
          zaxis = list(title = "Persistence"))
      )
  out.plot
  }

}


# if (plot) {
