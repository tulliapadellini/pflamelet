# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Persistence Flamelet
#'
#' Computes the Persistence Flamelet from a list of Persistence Diagrams. If the input is a list of data points observed at different scales, at each resolution, Persistence Diagrams are built
#' for the sub/superlevel set of
#' an arbitrary function computed on X, and then used to compute the Flamelet.
#'
#'@param X    a list of m persistence diagrams representing the scale-space family at different resolutions,
#'            or a list of m matrix/data.frame containing the pointcloud at different resolutions.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the \code{diag.fun}.
#'@param base.type a string specifying whether the Flamelet is built from Persistence Landscapes ("landscape") or Persistence Silhouettes ("silhouettes").
#'@param base.param the order k of the Flamelet (if \code{base.type=="landscape"}) or the power p of the Flamelet (if \code{base.type=="silhouette"}).
#'@param dimension the topological dimension of the flamelet (0 for connected components, 1 for loops, ...).
#'@param tseq a vector of dimension k containing the values at which the Flamelet function is evaluated for a fixed scale level.
#'@param diag.fun the function whose sub/super-level set define the persistent homology groups. Corresponds to the argument \code{FUN} of the \code{gridDiag} function in the TDA package.
#'@param sublevel a logical indicating whether the Persistent Homology should be computed on sub or superlevel set of the function given as \code{diag.fun}.
#'@param h.grid vector of dimension m containing bandwidths for the KDE, representing the scale parameter of the Flamelet.
#'@param lim 2-by-d matrix, where the i-th column contains the range of the grid over which the function specified in \code{diag.fun} is computed for the i-th variable.
#'@param by a scalar (or a vector if different values are selected for each dimension).
#'             specifying spaces between elements on the grid whose outernmost element are defined by \code{lim}.
#'@param scale a logical indicating whether or not the Persistence Diagrams have to be scaled to be in the same range (needed only for visualization purposes).
#'@param precomputed.diagram a logical indicating whether the Persistence Diagrams have to be computed or are given as input in the form of a list of diagrams.
#'@param info.message a logical denoting if progress messages should be printed.
#'@param band a logical indicating whether or not each Persistence Diagram should be cleaned by means of Bootstrap Bands
#'            (only possible when the Persistence Diagrams are computed within the function and are not given as an input).
#'@param B number of bootstrap repetitions needed to compute the confidence band over Persistence Diagrams.
#'@param alpha the confidence level of the bootstrap confidence bands.
#'@return a k-by-m matrix containing the Persistence Flamelet.
#'@references T. Padellini and P. Brutti (2017) Persistence Flamelets: multiscale Persistent Homology for kernel density exploration \url{https://arxiv.org/abs/1709.07097}
#'@examples
#'\donttest{
#'library(TDA)
#'xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'lim = cbind(Xlim, Ylim)
#'foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#'base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
#'                             tseq = seq(0, .75, length.out = 500))
#'}
#'@export
build.flamelet = function(X, base.type = "landscape", base.param = 1, dimension=1,
                           tseq, diag.fun = distFct, sublevel = TRUE, h.grid =NULL, lim = NULL, by=NULL,
                           scale = TRUE, precomputed.diagram = FALSE, info.message = FALSE, band = FALSE, B = 10, alpha = 0.95){

  if(!is.list(X)){
    f1 = function(h) {
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
    diagList = pbapply::pblapply(h.grid, f1)
    if(info.message) cat(" and now the Persistence Landscapes")

  } else {
    if(precomputed.diagram) diagList = X
    else{

      f2 = function(x) {
        diagr = gridDiag(x, diag.fun, lim = lim, by = by, sublevel = sublevel,
                         location = FALSE, printProgress = FALSE, maxdimension = dimension)$diagram
        return(diagr)
      }

      if(info.message) cat("Computing the Persistence Diagrams - hold on, this is the longest step...")
      diagList = pbapply::pblapply(X, f2)
      if(info.message) cat(" and now the Persistence Landscapes")

    }
  }

  if(base.type == "landscape") flamelet = sapply(diagList, landscape, dimension = dimension, KK = base.param, tseq = tseq)
  if(base.type == "silhouette") flamelet = sapply(diagList, silhouette, dimension = dimension, p = base.param, tseq = tseq)
  return(flamelet)

}


#' Bootstrap Band for Persistence Flamelet
#'
#' Computes a bootstrap band around the 0 of the Persistence Flamelet. If the input is a list of data points observed at different scales, at each resolution, Persistence Diagrams are built
#' for the sub/superlevel set of
#' an arbitrary function computed on X, and then used to compute the Flamelet.
#'
#'@param X    a list of m matrix/data.frame containing the pointcloud at different resolutions.
#'            If \code{X} is a n-by-d matrix or a data.frame containing a d-dimensional pointcloud,
#'            this function computes the Flamelet on the corresponding \code{diag.fun}.
#'@param B number of bootstrap repetitions needed to compute the confidence band over Persistence Diagrams.
#'@param alpha the confidence level of the bootstrap confidence bands.
#'@param base.type a string specifying whether the Flamelet is built from Persistence Landscapes ("landscape") or Persistence Silhouettes ("silhouettes").
#'@param base.param the order k of the Flamelet (if base.type=="landscape") or the power p of the Flamelet (if base.type=="silhouette").
#'@param dimension the topological dimension of the flamelet (0 for connected components, 1 for loops, ...)
#'@param tseq a vector of values at which the Flamelet function is evaluated for a fixed scale level.
#'@param diag.fun the function whose sub/super-level set define the persistent homology groups. Corresponds to the argument \code{FUN} of the \code{gridDiag} function in the TDA package.
#'@param sublevel a logical indicating whether the Persistent Homology should be computed on sub or superlevel set of the function given as \code{diag.fun}.
#'@param h.grid vector of bandwidths for the KDE, representing the scale parameter of the Flamelet.
#'@param lim 2-by-d matrix, where the i-th column contains the range of the grid over which the function specified in \code{diag.fun} is computed for the i-th variable.
#'@param by a scalar (or a vector if different values are selected for each dimension).
#'             specifying spaces between elements on the grid whose outernmost element are defined by \code{lim}.
#'@return The quantile of level alpha necessary to build the confidence band. More details can be found in Padellini (2017).
#'@references T. Padellini and P. Brutti (2017) Persistence Flamelets: multiscale Persistent Homology for kernel density exploration \url{https://arxiv.org/abs/1709.07097}
#'@examples
#'\donttest{
#'library(TDA)
#'xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'lim = cbind(Xlim, Ylim)
#'foo.band = flamelet.band(X = xx, B = 10, alpha = 0.95,
#'                   tseq = seq(0, .75, length.out = 500), diag.fun = kde,
#'                   h.grid = seq(0.01, 1, length.out = 40), lim = lim, by = by)
#'}
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
#' Plot the Persistence Flamelet in its original 3-dimensional form, or in a 2-dimensional projection.
#'
#'@param flamelet a kxm matrix corresponding to the Persistence Flamelet.
#'@param tseq a vector of length k containing of values at which the Flamelet function is evaluated for a fixed scale level.
#'@param scale.param a vector of lenght m corresponding to the values of the scale parameter at which the Flamelet has been evaluated.
#'@param flat a logical denoting whether the plot should be a 2-d projection of the Flamelet (\code{TRUE}) or a 3-d object (\code{FALSE}).
#'@param scale.name name of the scale parameter.
#'@param band scalar representing the confidence band for Persistence Flamelet. Only available when \code{flat = FALSE}.
#'@references T. Padellini and P. Brutti (2017) Persistence Flamelets: multiscale Persistent Homology for kernel density exploration \url{https://arxiv.org/abs/1709.07097}
#'@examples
#'\donttest{
#'library(TDA)
#'xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'lim = cbind(Xlim, Ylim)
#'foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#'base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
#'                             tseq = seq(0, .75, length.out = 500))
#'flamelet.plot(foo.flamelet, scale.param = seq(0.01, 1, length.out = 40),
#'              tseq = seq(0, .75, length.out = 500) )
#'
#'
#'foo.band = flamelet.band(X = xx, B = 10, alpha = 0.05,
#'                    tseq = seq(0, .75, length.out = 500), diag.fun = kde,
#'                    h.grid = seq(0.01, 1, length.out = 40), lim = lim, by = by)
#'
#'flamelet.plot(foo.flamelet, scale.param = seq(0.01, 1, length.out = 40),
#'              tseq = seq(0, .75, length.out = 500), band = foo.band)
#'}
#'@export
flamelet.plot = function(flamelet, scale.param, tseq, flat=FALSE,  scale.name ="Bandwidth", band = NULL){


  if(flat){
  image(y = scale.param , x = tseq, z= flamelet,
        col = viridis(100, option = "B"), ylab = scale.name, xlab = "(Birth + Death)/2")

  }else{
    out.plot = plot_ly(x = scale.param, y = tseq,  z = flamelet, colors = viridis(100, option = "B", end = .87)) %>%
      add_surface() %>%
      layout(
        scene = list(
          xaxis = list(title = scale.name),
          yaxis = list(title = "(Birth + Death)/2"),
          zaxis = list(title = "Persistence"))
      )
    if(!is.null(band)){
      upper.bound = flamelet + band
      lower.bound = flamelet - band
      lower.bound[lower.bound<0] <- 0
      out.plot = out.plot %>% add_surface(z = ~upper.bound, opacity = 0.5, showscale = FALSE) %>%
        add_surface(z = ~lower.bound, opacity = 0.5, showscale = FALSE)
    }

    out.plot
  }

}


#' Maximum Persistence Bandwidth
#'
#' Select the bandwidth that maximises the persistence of the scale-family, highlighting the most prominent topological feature of the KDE.
#'
#'@param flamelet a kxm matrix corresponding to the Persistence Flamelet.
#'@param scale.param a vector of lenght m corresponding to the values of the scale parameter at which the Flamelet \code{X} has been evaluated.
#'@return the value of the bandwidth corresponding to the most persistent feature of the Flamelet.
#'@references T. Padellini and P. Brutti (2017) Persistence Flamelets: multiscale Persistent Homology for kernel density exploration \url{https://arxiv.org/abs/1709.07097}
#'@examples
#'\donttest{
#'library(TDA)
#'xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'lim = cbind(Xlim, Ylim)
#'foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#'base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
#'                             tseq = seq(0, .75, length.out = 500))
#'mpband(flamelet = foo.flamelet, scale.param = seq(0.01, 1, length.out = 40))
#'}
#'@export
mpbandwidth = function(flamelet, scale.param){
  cc = apply(flamelet, 2, max)
  return(scale.param(which.max(cc)))

}



#' Clean Persistence Flamelet
#'
#' Remove all the values of the Persistence Flamelet which are not significantly different than 0
#'
#'@param flamelet a kxm matrix corresponding to the Persistence Flamelet.
#'@param band a scalar representing the confidence band.
#'@return the value of the bandwidth corresponding to the most persistent feature of the Flamelet.
#'@references T. Padellini and P. Brutti (2017) Persistence Flamelets: multiscale Persistent Homology for kernel density exploration \url{https://arxiv.org/abs/1709.07097}
#'@examples
#'\donttest{
#'library(TDA)
#'xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
#'Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
#'lim = cbind(Xlim, Ylim)
#'foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#'base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
#'                             tseq = seq(0, .75, length.out = 500))
#'foo.band = flamelet.band(X = xx, B = 10, alpha = 0.95,
#'                    tseq = seq(0, .75, length.out = 500), diag.fun = kde,
#'                    h.grid = seq(0.01, 1, length.out = 40), lim = lim, by = by)
#'new.flamelet = flamelet.clean(foo.flamelet, foo.band)
#'
#'}
#'@export
flamelet.clean = function(flamelet, band){
  idx = flamelet < band
  flamelet[idx] <- 0
  return(flamelet)
}


#' Two-sample permutation test
#'
#' Performs a permutation test to assess whether two samples of Persistence Flamelets
#' are likely be random draw from the same distribution or if they come from different
#' generating mechanisms, with p-value computed by means of bootstrap.
#'@param sample1 a k x m x n1 array corresponding to the Persistence Flamelet for the n1 subjects in sample 1
#'@param sample2 a k x m x n2 array corresponding to the Persistence Flamelet for the n2 subjects in sample 2
#'@param n.rep an integer representing the number of bootstrap iterations. If \code{NULL} only the test statistics
#'on the observed samples is returned
#'@param seed an integer specifying a seed for the random shuffling
#'@return a bootstrapped p-value
#'@references T. Padellini and P. Brutti (2017) Persistence Flamelets: multiscale Persistent Homology for kernel density exploration \url{https://arxiv.org/abs/1709.07097}
#'@examples
#'\donttest{
#'library(eegkit)
#'library(dplyr)
#'
#'# import data from the eegkit package
#' data("eegdata")   # electroencephalography data
#' data("eegcoord")  # location of the electrodes
#' # add eeg channel name as variable and select only the 2d projection of the electrode location
#' eegcoord <- mutate(eegcoord, channel = rownames(eegcoord)) %>% select(channel, xproj, yproj)
#'
#' # as EEG recordings are extremely unstable, they are typically averaged across repetitions
#' # here we average them across the 5 trials from each subject
#' eegmean <- eegdata %>% group_by(channel, time, subject) %>% summarise(mean = mean(voltage))
#' dim(eegmean) # 64 channels x 256 time points x 20 subjects
#'
#' subjects <- unique(eegdata$subject)
#' # subjects 1:10 are alcoholic, 11:20 are control
#'
#' # eegmean2 <- tapply(eegdata$voltage, list(eegdata$channel, eegdata$time, eegdata$subject), mean)
#' # Start by computing the list of Persistence Diagrams needed to build the flamelet for each subject
#' diag.list <- list()
#' t0 <- Sys.time()
#' for (sbj in subjects){
#'
#'   # select signal for one subject and then remove channels for which there are no coordinates
#'   sbj.data = dplyr::filter(eegmean, subject == sbj, !(channel %in% c("X", "Y", "nd")  ))
#'
#'   # add 2d projection of electrodes location
#'   sbj.data = left_join(sbj.data, eegcoord, by = "channel")
#'
#'   # scale data
#'   sbj.data[, c(4:6)] = scale(sbj.data[,c(4:6)])
#'
#'
#'   # dsucc.List = list()
#'
#'   diag.list.sbj = lapply(unique(sbj.data$time), function(time.idx){
#'     time.idx.data = filter(sbj.data, time == time.idx) %>% ungroup %>%
#'                     select(mean, xproj, yproj)
#'     time.idx.diag = ripsDiag(time.idx.data, maxdimension = 1, maxscale = 5,
#'                              library = "GUDHI", printProgress = F)
#'     return(time.idx.diag$diagram)
#'   })
#'
#'   diag.list[[which(sbj== subjects)]] = diag.list.sbj
#'
#'   print(paste("subject ", which(sbj == subjects), " of 20"))
#'
#'
#' }
#' t1 <- Sys.time()-t0
#' t1 # will take less than 5 minutes
#'
#' tseq <- seq(0, 5, length = 500) # consider 5 as it is the
#' # same value as maxscale (hence largest possible persistence)
#'
#' p_silh0 <- sapply(diag.list, FUN = build.flamelet,
#' base.type = "silhouette", dimension = 0,
#' tseq = tseq, precomputed.diagram = TRUE, simplify = 'array')
#'
#' prova = permutation_test(p_sih0[,,1:10], p_silh0[,,11:20],  n.rep = 10, seed = 1)
#'}
#'@export
permutation.test <- function(sample1, sample2,  n.rep = NULL, seed = NULL){
  # sample_1, sample_2 two samples of Persistence Flamelets

  n1 = dim(sample1)[3]
  n2 = dim(sample2)[3]
  nn = dim(sample1)[3] + dim(sample2)[3]

  stopifnot("Samples of Flamelets must have the same dimensions!" = !all(dim(sample1)[1:2] == dim(sample2)[1:2] ))
  sample_tot = abind::abind(sample1, sample2)


  flam1 = apply(sample_tot[,,1:n1],  c(1,2), mean)
  flam2 = apply(sample_tot[,,(n1+1):nn], c(1,2), mean)

  obs_diff = sum(abs(flam1-flam2))


  if(!is.null(n.rep)){
    if(!is.null(seed)) set.seed(seed)

    foo.fun = function(){
      shuffle_idx = sample(1:nn)
      sample_shuff = sample_tot[,,shuffle_idx]
      flam1 = apply(sample_shuff[,,1:n1],   c(1,2),mean)
      flam2 = apply(sample_shuff[,,(n1+1):nn], c(1,2),mean)

      (sum(abs(flam1-flam2)))}

    shuff_diff = pbapply::pbreplicate(n.rep, expr = foo.fun())

    return(mean(shuff_diff>obs_diff))

  } else {
    return(obs_diff)
  }
}

