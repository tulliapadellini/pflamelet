# library(TDA)
# xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
# Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
# lim = cbind(Xlim, Ylim)
#
# build.flamelet = function(X, base.type = "landscape", base.param = 1, dimension=1,
#                           tseq, diag.fun = distFct, h.grid =NULL, lim = NULL, by=NULL,
#                           scale = TRUE, precomputed.diagram = FALSE, sublevel = TRUE)
#
#
# foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
#                               base.type = "landscape", dimension = 1,base.param = 1,
#                               lim = lim, by = by,
#                               tseq = seq(0, .12, length.out = 500), scale = F)
#
# image(foo.flamelet)
#
#
# flamelet.band = function(X, B, alpha, base.type = "landscape", base.param = 1, dimension=1,
#                          tseq, diag.fun = distFct, h.grid =NULL, lim = NULL, by=NULL,
#                          sublevel = TRUE)
#
#
#
# kk = flamelet.band(xx, 10, 0.05, h.grid = seq(0.01, 1, length.out = 40),
#                    base.type = "landscape", dimension = 1,base.param = 2,
#                    lim = lim, by = by,
#                    tseq = seq(0, .12, length.out = 500))
#
#
#
# kk2 = flamelet.band(xx, 100, 0.05, h.grid = seq(0.01, 1, length.out = 40),
#                    base.type = "landscape", dimension = 1,base.param = 1,
#                    lim = lim, by = by,
#                    tseq = seq(0, .12, length.out = 500))
#
#
# kk3 = flamelet.band(xx, 100, 0.95, h.grid = seq(0.01, 1, length.out = 40),
#                     base.type = "landscape", dimension = 1,base.param = 1,
#                     lim = lim, by = by,
#                     tseq = seq(0, .12, length.out = 500))
#
#
# kk
# kk2
# sum((foo.flamelet - kk)>0)
# image((foo.flamelet - kk)>0)
#
# sum((foo.flamelet - kk3)>0)
# image((foo.flamelet - kk2)>0)
#
#
# image(foo.flamelet)
# plot_ly(z = foo.flamelet, type = "surface")
