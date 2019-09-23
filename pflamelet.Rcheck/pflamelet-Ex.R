pkgname <- "pflamelet"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "pflamelet-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('pflamelet')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("build.flamelet")
### * build.flamelet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: build.flamelet
### Title: Persistence Flamelet Function
### Aliases: build.flamelet

### ** Examples


## library(TDA)
## xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
## Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
## lim = cbind(Xlim, Ylim)
## foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
## base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("build.flamelet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("flamelet.plot")
### * flamelet.plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: flamelet.plot
### Title: Plot Persistence Flamelet
### Aliases: flamelet.plot

### ** Examples


## library(TDA)
## xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
## Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
## lim = cbind(Xlim, Ylim)
## foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
## base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("flamelet.plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
