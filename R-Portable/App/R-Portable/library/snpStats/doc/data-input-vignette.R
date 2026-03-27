### R code from vignette source 'data-input-vignette.Rnw'

###################################################
### code chunk number 1: lib
###################################################
require(snpStats)


###################################################
### code chunk number 2: sysfile
###################################################
pedfile <- system.file("extdata/sample.ped.gz", package="snpStats")
pedfile
infofile <- system.file("extdata/sample.info", package="snpStats")


###################################################
### code chunk number 3: redpedfile
###################################################
sample <- read.pedfile(pedfile, snps=infofile)


###################################################
### code chunk number 4: sample1
###################################################
sample$genotypes
col.summary(sample$genotypes)$MAF
head(col.summary(sample$genotypes))


###################################################
### code chunk number 5: sample2
###################################################
head(sample$fam)


###################################################
### code chunk number 6: sample3
###################################################
head(sample$map)


###################################################
### code chunk number 7: plink
###################################################
fam <- system.file("extdata/sample.fam", package="snpStats")
bim <- system.file("extdata/sample.bim", package="snpStats")
bed <- system.file("extdata/sample.bed", package="snpStats")
sample <- read.plink(bed, bim, fam)


###################################################
### code chunk number 8: plinkout
###################################################
sample$genotypes
col.summary(sample$genotypes)$MAF
head(sample$fam)
head(sample$map)


###################################################
### code chunk number 9: plinkselect
###################################################
subset <- read.plink(bed, bim, fam, select.snps=6:10)
subset$genotypes
col.summary(subset$genotypes)$MAF
subset$map


###################################################
### code chunk number 10: longfile
###################################################
longfile <- system.file("extdata/sample-long.gz", package="snpStats")
longfile


###################################################
### code chunk number 11: longlist
###################################################
cat(readLines(longfile, 5), sep="\n")


###################################################
### code chunk number 12: readlong
###################################################
gdata <- read.long(longfile, 
   fields=c(snp=1, sample=2, genotype=3, confidence=4),
   gcodes=c("1", "2", "3"), 
   threshold=0.95)
gdata
summary(gdata)


###################################################
### code chunk number 13: readlongallele
###################################################
allelesfile <- system.file("extdata/sample-long-alleles.gz", package="snpStats")
cat(readLines(allelesfile, 5), sep="\n")
gdata <- read.long(allelesfile, 
   fields=c(snp=1, sample=2, allele.A=3, allele.B=4, confidence=5),
   threshold=0.95)
gdata
gdata$genotypes
gdata$alleles


