### R code from vignette source 'pca-vignette.Rnw'

###################################################
### code chunk number 1: get-data
###################################################
require(snpStats)
data(for.exercise)
controls <- rownames(subject.support)[subject.support$cc==0]
use <- seq(1, ncol(snps.10), 10)
ctl.10 <- snps.10[controls,use]


###################################################
### code chunk number 2: xxt-matrix
###################################################
xxmat <- xxt(ctl.10, correct.for.missing=FALSE)


###################################################
### code chunk number 3: eigen
###################################################
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]
evals <- evv$values[1:5]
evals


###################################################
### code chunk number 4: pc-one
###################################################
pop <- subject.support[controls,"stratum"]
par(mfrow=c(1,2))
boxplot(pcs[,1]~pop)
boxplot(pcs[,2]~pop)


###################################################
### code chunk number 5: pre-multiply
###################################################
btr <- snp.pre.multiply(ctl.10, diag(1/sqrt(evals)) %*% t(pcs))


###################################################
### code chunk number 6: post-multiply
###################################################
pcs <- snp.post.multiply(snps.10[,use], t(btr))


###################################################
### code chunk number 7: testing
###################################################
cc <- subject.support$cc
uncorrected <- single.snp.tests(cc, snp.data=snps.10)
corrected <- snp.rhs.tests(cc~pcs[,1], snp.data=snps.10)
par(mfrow=c(1,2),cex.sub=0.85)
qq.chisq(chi.squared(uncorrected,1), df=1)
qq.chisq(chi.squared(corrected), df=1)


