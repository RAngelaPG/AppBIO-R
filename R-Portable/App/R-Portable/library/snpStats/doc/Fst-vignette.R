### R code from vignette source 'Fst-vignette.Rnw'

###################################################
### code chunk number 1: lib
###################################################
require(snpStats)


###################################################
### code chunk number 2: pairwise
###################################################
data(for.exercise)
f1 <- Fst(snps.10, subject.support$stratum, pairwise=TRUE)
weighted.mean(f1$Fst, f1$weight)


###################################################
### code chunk number 3: aov
###################################################
f2 <- Fst(snps.10, subject.support$stratum, pairwise=FALSE)
weighted.mean(f2$Fst, f2$weight)


###################################################
### code chunk number 4: groupsize
###################################################
table(subject.support$stratum)


