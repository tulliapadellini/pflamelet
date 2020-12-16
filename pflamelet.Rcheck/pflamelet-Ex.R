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
### Title: Persistence Flamelet
### Aliases: build.flamelet

### ** Examples

## No test: 
library(TDA)
xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
lim = cbind(Xlim, Ylim)
foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
                            tseq = seq(0, .75, length.out = 500))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("build.flamelet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("flamelet.band")
### * flamelet.band

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: flamelet.band
### Title: Bootstrap Band for Persistence Flamelet
### Aliases: flamelet.band

### ** Examples

## No test: 
library(TDA)
xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
lim = cbind(Xlim, Ylim)
foo.band = flamelet.band(X = xx, B = 10, alpha = 0.95,
                  tseq = seq(0, .75, length.out = 500), diag.fun = kde,
                  h.grid = seq(0.01, 1, length.out = 40), lim = lim, by = by)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("flamelet.band", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("flamelet.clean")
### * flamelet.clean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: flamelet.clean
### Title: Clean Persistence Flamelet
### Aliases: flamelet.clean

### ** Examples

## No test: 
library(TDA)
xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
lim = cbind(Xlim, Ylim)
foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
                            tseq = seq(0, .75, length.out = 500))
foo.band = flamelet.band(X = xx, B = 10, alpha = 0.95,
                   tseq = seq(0, .75, length.out = 500), diag.fun = kde,
                   h.grid = seq(0.01, 1, length.out = 40), lim = lim, by = by)
new.flamelet = flamelet.clean(foo.flamelet, foo.band)

## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("flamelet.clean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("flamelet.plot")
### * flamelet.plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: flamelet.plot
### Title: Plot Persistence Flamelet
### Aliases: flamelet.plot

### ** Examples

## No test: 
library(TDA)
xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
lim = cbind(Xlim, Ylim)
foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
                            tseq = seq(0, .75, length.out = 500))
flamelet.plot(foo.flamelet, scale.param = seq(0.01, 1, length.out = 40),
             tseq = seq(0, .75, length.out = 500) )


foo.band = flamelet.band(X = xx, B = 10, alpha = 0.05,
                   tseq = seq(0, .75, length.out = 500), diag.fun = kde,
                   h.grid = seq(0.01, 1, length.out = 40), lim = lim, by = by)

flamelet.plot(foo.flamelet, scale.param = seq(0.01, 1, length.out = 40),
             tseq = seq(0, .75, length.out = 500), band = foo.band)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("flamelet.plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mpbandwidth")
### * mpbandwidth

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mpbandwidth
### Title: Maximum Persistence Bandwidth
### Aliases: mpbandwidth

### ** Examples

## No test: 
library(TDA)
xx = rbind(circleUnif(50, 1), circleUnif(50, 1.5) + 3)
Xlim = c(-1, 5);  Ylim = c(-1, 5);  by = 0.05
lim = cbind(Xlim, Ylim)
foo.flamelet = build.flamelet(X = xx, h.grid = seq(0.01, 1, length.out = 40),
base.type = "landscape", dimension = 1,base.param = 1, lim = lim, by = by,
                            tseq = seq(0, .75, length.out = 500))
mpband(flamelet = foo.flamelet, scale.param = seq(0.01, 1, length.out = 40))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mpbandwidth", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("permutation.test")
### * permutation.test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: permutation.test
### Title: Two-sample permutation test
### Aliases: permutation.test

### ** Examples

## No test: 
library(eegkit)
library(dplyr)

# import data from the eegkit package
data("eegdata")   # electroencephalography data
data("eegcoord")  # location of the electrodes
# add eeg channel name as variable and select only the 2d projection of the electrode location
eegcoord <- mutate(eegcoord, channel = rownames(eegcoord)) %>% select(channel, xproj, yproj)

# as EEG recordings are extremely unstable, they are typically averaged across repetitions
# here we average them across the 5 trials from each subject
eegmean <- eegdata %>% group_by(channel, time, subject) %>% summarise(mean = mean(voltage))
dim(eegmean) # 64 channels x 256 time points x 20 subjects

subjects <- unique(eegdata$subject)
# subjects 1:10 are alcoholic, 11:20 are control

# eegmean2 <- tapply(eegdata$voltage, list(eegdata$channel, eegdata$time, eegdata$subject), mean)
# Start by computing the list of Persistence Diagrams needed to build the flamelet for each subject
diag.list <- list()
t0 <- Sys.time()
for (sbj in subjects){

  # select signal for one subject and then remove channels for which there are no coordinates
  sbj.data = dplyr::filter(eegmean, subject == sbj, !(channel %in% c("X", "Y", "nd")  ))

  # add 2d projection of electrodes location
  sbj.data = left_join(sbj.data, eegcoord, by = "channel")

  # scale data
  sbj.data[, c(4:6)] = scale(sbj.data[,c(4:6)])


  # dsucc.List = list()

  diag.list.sbj = lapply(unique(sbj.data$time), function(time.idx){
    time.idx.data = filter(sbj.data, time == time.idx) %>% ungroup %>%
                    select(mean, xproj, yproj)
    time.idx.diag = ripsDiag(time.idx.data, maxdimension = 1, maxscale = 5,
                             library = "GUDHI", printProgress = F)
    return(time.idx.diag$diagram)
  })

  diag.list[[which(sbj== subjects)]] = diag.list.sbj

  print(paste("subject ", which(sbj == subjects), " of 20"))


}
t1 <- Sys.time()-t0
t1 # will take less than 5 minutes

tseq <- seq(0, 5, length = 500) # consider 5 as it is the
# same value as maxscale (hence largest possible persistence)

p_silh0 <- sapply(diag.list, FUN = build.flamelet,
base.type = "silhouette", dimension = 0,
tseq = tseq, precomputed.diagram = TRUE, simplify = 'array')

prova = permutation_test(p_sih0[,,1:10], p_silh0[,,11:20],  n.rep = 10, seed = 1)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("permutation.test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
