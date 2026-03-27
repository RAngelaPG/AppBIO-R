### R code from vignette source 'residuetheorem.Rnw'

###################################################
### code chunk number 1: requirepackage
###################################################
require(elliptic,quietly=TRUE)


###################################################
### code chunk number 2: residuetheorem.Rnw:159-160
###################################################
integrate.segments(exp, c(0, 1, 1+1i, 1i), close=TRUE)


###################################################
### code chunk number 3: residuetheorem.Rnw:173-176
###################################################
analytic <- exp(1)*(exp(1i)-1)
numeric  <- integrate.segments(exp, c(1, 1+1i), close=FALSE)
c(analytic=analytic, numeric=numeric, difference=analytic-numeric)


###################################################
### code chunk number 4: residuetheorem.Rnw:202-208
###################################################
u     <- function(x){exp(pi*2i*x)}
udash <- function(x){pi*2i * exp(pi*2i*x)}

analytic <- pi*2i
numeric  <- integrate.contour(function(z){1/z}, u, udash)
c(analytic=analytic, numeric=numeric, difference=analytic-numeric)


###################################################
### code chunk number 5: showhypergeofail
###################################################
library("hypergeo")
z0 <- 1/2 + sqrt(3)/2i
f <- function(z){hypergeo_powerseries(1/2, 1/3, 1/5, z)}
f(z0)


###################################################
### code chunk number 6: residuetheorem.Rnw:272-276
###################################################
r <- 0.1 # radius of contour
u <- function(x){z0 + r*exp(pi * 2i * x)}
udash <- function(x){r * pi * (0+2i) * exp(pi * 2i * x)}
(val_residue <- integrate.contour(function(z){f(z) / (z-z0)}, u, udash) / (pi*2i))


###################################################
### code chunk number 7: residuetheorem.Rnw:283-285
###################################################
(val_gosper <- hypergeo_gosper(1/2, 1/3, 1/5, z0))
abs(val_gosper - val_residue)


