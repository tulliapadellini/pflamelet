xx = rbind(circleUnif(100, 1), circleUnif(100, 1.5) + 3)
plot(xx)

Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
lim = cbind(Xlim, Ylim)


foo.diag  = gridDiag(xx, kde, lim = lim, by = by, sublevel = FALSE, h = 0.3)

plot(foo.diag$diagram)

build.flamelet2 = function(X, base.type = "landscape", base.param = 1, dimension=1,
                          tseq, diag.fun = distFct, h.grid =NULL, lim = NULL, by=NULL,
                          scale = TRUE, precomputed.diagram = FALSE, sublevel = TRUE, info.message = FALSE, band = FALSE, B = 10, alpha = 0.95){

  if(!is.list(X)){
    f = function(h) {
      diagr = gridDiag(X, kde, lim = lim, by = by, sublevel = FALSE,
                       location = FALSE, printProgress = FALSE, h = h, maxdimension = dimension)$diagram

      if(band) {
        cc <- bootstrapDiagram(X, kde, lim = lim, by = by, sublevel = FALSE, B = B,
                             alpha = 1-alpha/length(h.grid), dimension = dimension, printProgress = FALSE, h = h)


        idx = apply(diagr[,2:3], 1, diff)>2*cc
        diagr = diagr[idx,]

      }

      if(scale && nrow(diagr)!=0) diagr[,2:3] = diagr[,2:3]/diagr[1,3]
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




foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 0.7, length.out = 20),
                              base.type = "landscape", dimension = 1,base.param = 1,
                              lim = lim, by = by,
                              tseq = seq(0, 1, length.out = 500), scale = TRUE, band = TRUE)


image(foo.flamelet)
plot_ly(z = foo.flamelet, type = "surface")
