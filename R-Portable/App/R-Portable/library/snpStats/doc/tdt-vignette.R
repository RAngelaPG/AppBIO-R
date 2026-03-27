### R code from vignette source 'tdt-vignette.Rnw'

###################################################
### code chunk number 1: family-data
###################################################
require(snpStats)
data(families)
genotypes
head(pedData)


###################################################
### code chunk number 2: mis-inheritances
###################################################
mis <- misinherits(data=pedData, snp.data=genotypes)
dim(mis)


###################################################
### code chunk number 3: per-subj-snp
###################################################
per.subj <- apply(mis, 1, sum, na.rm=TRUE)
per.snp <- apply(mis, 2, sum, na.rm=TRUE)
par(mfrow = c(1, 2))
hist(per.subj,main='Histogram per Subject', xlab='Subject')
hist(per.snp,main='Histogram per SNP', xlab='SNP')


###################################################
### code chunk number 4: per-family
###################################################
fam <- pedData[rownames(mis), "familyid"]
per.fam <- tapply(per.subj, fam, sum)
par(mfrow = c(1, 1))
hist(per.fam, main='Histogram per Family', xlab='Family')


###################################################
### code chunk number 5: tdt-tests
###################################################
tests <- tdt.snp(data = pedData, snp.data = genotypes)
cbind(p.values.1df = p.value(tests, 1),
      p.values.2df = p.value(tests, 2))
qq.chisq(chi.squared(tests, 1), df = 1)


