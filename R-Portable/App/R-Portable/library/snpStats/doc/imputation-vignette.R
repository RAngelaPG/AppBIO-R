### R code from vignette source 'imputation-vignette.Rnw'

###################################################
### code chunk number 1: init
###################################################
library(snpStats)
library(hexbin)
data(for.exercise)


###################################################
### code chunk number 2: select
###################################################
training <- sample(1000, 200)
select <- seq(1, ncol(snps.10),2) 
missing <- snps.10[training, select]
present <- snps.10[training, -select] 
missing 
present 


###################################################
### code chunk number 3: target
###################################################
target <- snps.10[-training, -select]
target


###################################################
### code chunk number 4: imputation-vignette.Rnw:91-93
###################################################
lost <- snps.10[-training, select]
lost


###################################################
### code chunk number 5: positions
###################################################
pos.miss <- snp.support$position[select]
pos.pres <- snp.support$position[-select]


###################################################
### code chunk number 6: rules
###################################################
rules <- snp.imputation(present, missing, pos.pres, pos.miss)


###################################################
### code chunk number 7: rule1
###################################################
rules[1:10]


###################################################
### code chunk number 8: rule2
###################################################
rules[c('rs11253563', 'rs2379080')]


###################################################
### code chunk number 9: summary
###################################################
summary(rules)


###################################################
### code chunk number 10: ruleplot
###################################################
plot(rules)


###################################################
### code chunk number 11: imptest
###################################################
imp <- single.snp.tests(cc, stratum, data=subject.support,
                        snp.data=target, rules=rules)


###################################################
### code chunk number 12: realtest
###################################################
obs <- single.snp.tests(cc, stratum, data=subject.support, snp.data=lost)


###################################################
### code chunk number 13: compare
###################################################
logP.imp <- -log10(p.value(imp, df=1))
logP.obs <- -log10(p.value(obs, df=1))
hb <- hexbin(logP.obs, logP.imp, xbin=50)
sp <- plot(hb)
hexVP.abline(sp$plot.vp, 0, 1, col="black")


###################################################
### code chunk number 14: best
###################################################
use <- imputation.r2(rules)>0.9
hb <- hexbin(logP.obs[use], logP.imp[use], xbin=50)
sp <- plot(hb)
hexVP.abline(sp$plot.vp, 0, 1, col="black")


###################################################
### code chunk number 15: rsqmaf
###################################################
hb <- hexbin(imputation.maf(rules), imputation.r2(rules), xbin=50)
sp <- plot(hb)


###################################################
### code chunk number 16: imptest-rhs
###################################################
imp2 <- snp.rhs.tests(cc~strata(stratum), family="binomial", 
                      data=subject.support, snp.data=target, rules=rules)
logP.imp2 <- -log10(p.value(imp2))
hb <- hexbin(logP.obs, logP.imp2, xbin=50)
sp <- plot(hb)
hexVP.abline(sp$plot.vp, 0, 1, col="black")


###################################################
### code chunk number 17: impstore
###################################################
imputed <- impute.snps(rules, target, as.numeric=FALSE)


###################################################
### code chunk number 18: uncert1
###################################################
plotUncertainty(imputed[, "rs4880568"])


###################################################
### code chunk number 19: uncert2
###################################################
plotUncertainty(imputed[, "rs2050968"])


###################################################
### code chunk number 20: imptest2
###################################################
imp3 <- single.snp.tests(cc, stratum, data=subject.support,
                        snp.data=imputed, uncertain=TRUE)


###################################################
### code chunk number 21: imp3
###################################################
imp3[1:5]
imp[1:5]


###################################################
### code chunk number 22: mach
###################################################
path <- system.file("extdata/mach1.out.mlprob.gz", package="snpStats")
mach <- read.mach(path)
plotUncertainty(mach[,50])


###################################################
### code chunk number 23: class-imp-obs
###################################################
class(imp)


###################################################
### code chunk number 24: save-scores
###################################################
obs <- single.snp.tests(cc, stratum, data=subject.support, snp.data=missing, 
                        score=TRUE)
imp <- single.snp.tests(cc, stratum, data=subject.support,
                        snp.data=target, rules=rules, score=TRUE)


###################################################
### code chunk number 25: imputation-vignette.Rnw:345-347
###################################################
class(obs)
class(imp)


###################################################
### code chunk number 26: pool
###################################################
both <- pool(obs, imp)
class(both)
both[1:5]


###################################################
### code chunk number 27: pool-score
###################################################
both <- pool(obs, imp, score=TRUE)
class(both)


###################################################
### code chunk number 28: sign
###################################################
table(effect.sign(obs))


###################################################
### code chunk number 29: switch
###################################################
effect.sign(obs)[1:6]
sw.obs <- switch.alleles(obs, 1:3)
class(sw.obs)
effect.sign(sw.obs)[1:6]


