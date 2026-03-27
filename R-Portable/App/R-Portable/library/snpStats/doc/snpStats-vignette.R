### R code from vignette source 'snpStats-vignette.Rnw'

###################################################
### code chunk number 1: init
###################################################
require(snpStats)
require(hexbin)
data(for.exercise)


###################################################
### code chunk number 2: snpStats-vignette.Rnw:105-106
###################################################
show(snps.10)


###################################################
### code chunk number 3: snpStats-vignette.Rnw:112-113
###################################################
summary(snp.support)


###################################################
### code chunk number 4: snpStats-vignette.Rnw:124-125
###################################################
summary(subject.support)


###################################################
### code chunk number 5: snpStats-vignette.Rnw:146-147
###################################################
summary(snps.10)


###################################################
### code chunk number 6: snpStats-vignette.Rnw:152-154
###################################################
snpsum <- col.summary(snps.10)
summary(snpsum)


###################################################
### code chunk number 7: plot-snpsum
###################################################
par(mfrow = c(1, 2))
hist(snpsum$MAF)
hist(snpsum$z.HWE)


###################################################
### code chunk number 8: sample-qc
###################################################
sample.qc <- row.summary(snps.10)
summary(sample.qc)


###################################################
### code chunk number 9: plot-outliners-qc
###################################################
par(mfrow = c(1, 1))
plot(sample.qc)


###################################################
### code chunk number 10: outliers
###################################################
use <- sample.qc$Heterozygosity>0
snps.10 <- snps.10[use, ]
subject.support <- subject.support[use, ]


###################################################
### code chunk number 11: if-case-control
###################################################
if.case <- subject.support$cc == 1
if.control <- subject.support$cc == 0


###################################################
### code chunk number 12: sum-case-control
###################################################
sum.cases <- col.summary(snps.10[if.case, ])
sum.controls <- col.summary(snps.10[if.control, ])


###################################################
### code chunk number 13: plot-summaries
###################################################
hb <- hexbin(sum.controls$Call.rate, sum.cases$Call.rate, xbin=50)
sp <- plot(hb)
hexVP.abline(sp$plot.vp, 0, 1, col="black")


###################################################
### code chunk number 14: plot-freqs
###################################################
sp <- plot(hexbin(sum.controls$MAF, sum.cases$MAF, xbin=50))
hexVP.abline(sp$plot.vp, 0, 1, col="white")


###################################################
### code chunk number 15: tests
###################################################
tests <- single.snp.tests(cc, data = subject.support, snp.data = snps.10)


###################################################
### code chunk number 16: sum-tests
###################################################
summary(tests)


###################################################
### code chunk number 17: use
###################################################
use <- snpsum$MAF > 0.01 & snpsum$z.HWE^2 < 200


###################################################
### code chunk number 18: sum-use
###################################################
sum(use)


###################################################
### code chunk number 19: subset-tests
###################################################
tests <- tests[use]
position <- snp.support[use, "position"]


###################################################
### code chunk number 20: plot-tests
###################################################
p1 <- p.value(tests, df=1)
plot(hexbin(position, -log10(p1), xbin=50))


###################################################
### code chunk number 21: qqplot
###################################################
chi2 <- chi.squared(tests, df=1)
qq.chisq(chi2,  df = 1)


###################################################
### code chunk number 22: more-tests
###################################################
tests <- single.snp.tests(cc, stratum, data = subject.support,
     snp.data = snps.10)
tests <- tests[use]
p1 <- p.value(tests, df = 1)
plot(hexbin(position, -log10(p1), xbin=50))


###################################################
### code chunk number 23: more-tests-qq
###################################################
chi2 <- chi.squared(tests, df=1)
qq.chisq(chi2, df = 1)


###################################################
### code chunk number 24: ord
###################################################
ord <- order(p1)
top10 <- ord[1:10]
top10


###################################################
### code chunk number 25: top-10
###################################################
names <- tests@snp.names
p1[top10]
names[top10]
position[top10]


###################################################
### code chunk number 26: top10-local
###################################################
posord <- order(position)
position <- position[posord]
names <- names[posord]
local <- names[position > 9.6e+07 & position < 9.8e+07]


###################################################
### code chunk number 27: top1
###################################################
cc <- subject.support$cc
stratum <- subject.support$stratum
top <- as(snps.10[, "rs870041"], "numeric")
glm(cc ~ top + stratum, family = "binomial")


###################################################
### code chunk number 28: top2
###################################################
top2 <- as(snps.10[, "rs10882596"], "numeric")
fit <- glm(cc ~ top2 + stratum, family = "binomial")
summary(fit)


###################################################
### code chunk number 29: estimates
###################################################
localest <- snp.rhs.estimates(cc~stratum, family="binomial", sets=local,
                              data=subject.support, snp.data=snps.10)


###################################################
### code chunk number 30: list-estimates
###################################################
localest[1:5]
localest["rs10882596"]


###################################################
### code chunk number 31: fast-estimates
###################################################
allest <- snp.rhs.estimates(cc~strata(stratum), family="binomial", sets=use,
                              data=subject.support, snp.data=snps.10)
length(allest)


###################################################
### code chunk number 32: check-estimates
###################################################
allest["rs10882596"]


###################################################
### code chunk number 33: blocks
###################################################
blocks <- split(posord, cut(position, seq(100000, 135300000, 20000)))
bsize <- sapply(blocks, length)
table(bsize)
blocks <- blocks[bsize>0]  


###################################################
### code chunk number 34: twentyfive
###################################################
posord[1:20]
blocks[1:5]


###################################################
### code chunk number 35: blocks-use
###################################################
snps.use <- snps.10[, use]
remove(snps.10)


###################################################
### code chunk number 36: mtests
###################################################
mtests <- snp.rhs.tests(cc ~ stratum, family = "binomial", 
     data = subject.support, snp.data = snps.use, tests = blocks)
summary(mtests)


###################################################
### code chunk number 37: plot-mtests
###################################################
pm <- p.value(mtests)
plot(hexbin(-log10(pm), xbin=50))


###################################################
### code chunk number 38: qqplot-mtests
###################################################
qq.chisq(-2 * log(pm), df = 2)


