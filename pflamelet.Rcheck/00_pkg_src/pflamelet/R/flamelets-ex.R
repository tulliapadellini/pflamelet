library(TDA)

# these are for visualization
library(plotly)
library(viridis)


build.flamelet = function(X, k, tseq, dimension, progr = FALSE, plot = FALSE, h.grid =NULL, lim = NULL, by=NULL){

  ## arguments always needed:

  # X:         matrix,dataframe / list - the data, either as a pointcloud, either as a list
  #            of persistence diagrams
  # k:         integer - the order of the flamelet
  # dimension: integer - the dimension of the flamelet
  # progr:     logical - if TRUE
  # plot:      logical - if TRUE the function compute the 3d plot

  ## if X is a matrix or a dataframe the function computes the flamelets of the
  ## kernel density estimator computed on X, in this case, the following
  ## additional arguments have to be specified:

  # h.grid:    numeric - vector of bandwidths for the KDE
  # lim:       range of each dimension of the grid over which KDE is computed
  # by:        numeric - scalar or (vector if different values are selected for each dimension)
  #            specifying spaces between elements on the grid



  if(!is.list(X)){
    f = function(h) {
      diagr = gridDiag(X, kde, lim = lim, by = by, sublevel = FALSE,
                       location = FALSE, printProgress = progr, h = h, maxdimension = dimension)$diagram
      diagr[,2:3] = diagr[,2:3]/diagr[1,3]
      return(diagr)
    }

    cat("Computing the Persistence Diagrams - hold on, this is the longest step...")
    diagListh = lapply(h.grid, f)
    cat(" and now the Persistence Landscapes")

  } else {
    diagListh = X
  }

   landSurf = sapply(diagListh, landscape, dimension = dimension, KK = k, tseq = tseq)


   if (plot) {
     three.d = plot_ly(x = h.grid, y = tseq,  z = landSurf, colors = viridis(100, option = "B", end = .87)) %>%
       add_surface() %>%
       layout(
         scene = list(
           xaxis = list(title = "Bandwidth"),
           yaxis = list(title = "(Birth + Death)/2"),
           zaxis = list(title = "Persistence"))
       )

     return(list(three.d, landSurf) )
   } else   return(landSurf)

}


### 0. Create the data
#
# xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
# plot(xx, pch = 20)
#
# Xlim <- c(-1, 5);  Ylim <- c(-1, 5);  by <- 0.05
# lim = cbind(Xlim, Ylim)
#
# ### 1. Compute the Flamelets
#
# foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40), k = 1, lim = lim, by = by,
#                               tseq = seq(0, .75, length.out = 500), progr = F, plot = TRUE)
#
#
# ### 2. Visualize the Flamelets
#
#
# # 3d plot
# foo.flamelet[[1]]
#
# # 2d plot
# image(y = seq(0.01, 1, length.out = 40) , x = seq(0, .75, length.out = 500), z= foo.flamelet[[2]],
#       col = viridis(100, option = "B"), ylab = "bandwidth", xlab = "(birth + death)/2")
#
