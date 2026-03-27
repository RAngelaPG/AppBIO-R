### R code from vignette source 'ReadingFiles.Rnw'

###################################################
### code chunk number 1: ReadingFiles.Rnw:30-31
###################################################
options(width = 80, prompt = "> ")


###################################################
### code chunk number 2: ReadingFiles.Rnw:107-112
###################################################
library(pegas)
x <- read.loci("toto", header = FALSE)
x
print(x, details = TRUE)
class(x)


###################################################
### code chunk number 3: ReadingFiles.Rnw:130-133
###################################################
y <- read.loci("titi", header = FALSE, allele.sep = "-")
print(y, details = TRUE)
identical(x, y)


###################################################
### code chunk number 4: ReadingFiles.Rnw:138-139
###################################################
args(read.loci)


###################################################
### code chunk number 5: ReadingFiles.Rnw:166-167
###################################################
print(read.loci("tutu", FALSE), TRUE)


###################################################
### code chunk number 6: ReadingFiles.Rnw:183-186
###################################################
X <- read.loci("tyty")
print(X, TRUE)
summary(X)


###################################################
### code chunk number 7: ReadingFiles.Rnw:204-207
###################################################
z <- read.loci("tata", loci.sep = "\t", col.loci = 2:3, col.pop = 4, row.names = 1)
z
print(z, details = TRUE)


###################################################
### code chunk number 8: ReadingFiles.Rnw:213-214
###################################################
getAlleles(z)


###################################################
### code chunk number 9: ReadingFiles.Rnw:219-220
###################################################
attr(z, "locicol")


###################################################
### code chunk number 10: ReadingFiles.Rnw:225-226
###################################################
str(z)


###################################################
### code chunk number 11: ReadingFiles.Rnw:234-235
###################################################
args(read.vcf)


