### R code from vignette source 'ld-vignette.Rnw'

###################################################
### code chunk number 1: lib
###################################################
require(snpStats)
require(hexbin)


###################################################
### code chunk number 2: data
###################################################
data(ld.example)


###################################################
### code chunk number 3: showgt
###################################################
ceph.1mb
yri.1mb


###################################################
### code chunk number 4: showsp
###################################################
head(support.ld)


###################################################
### code chunk number 5: ldstats
###################################################
ld.ceph <- ld(ceph.1mb, stats=c("D.prime", "R.squared"), depth=100)
ld.yri <- ld(yri.1mb, stats=c("D.prime", "R.squared"), depth=100)


###################################################
### code chunk number 6: image1 (eval = FALSE)
###################################################
## image(ld.ceph$D.prime, lwd=0)


###################################################
### code chunk number 7: image1a
###################################################
print(image(ld.ceph$D.prime, lwd=0))


###################################################
### code chunk number 8: image2 (eval = FALSE)
###################################################
## image(ld.yri$D.prime, lwd=0)


###################################################
### code chunk number 9: image2a
###################################################
print(image(ld.yri$D.prime, lwd=0))


###################################################
### code chunk number 10: quartiles
###################################################
quantile(ld.ceph$D.prime@x, na.rm=TRUE)
quantile(ld.yri$D.prime@x, na.rm=TRUE)


###################################################
### code chunk number 11: colors
###################################################
spectrum <- rainbow(10, start=0, end=1/6)[10:1]


###################################################
### code chunk number 12: imagecol (eval = FALSE)
###################################################
## image(ld.ceph$D.prime, lwd=0, cuts=9, col.regions=spectrum, colorkey=TRUE)


###################################################
### code chunk number 13: imagecola
###################################################
print(image(ld.ceph$D.prime, lwd=0, cuts=9, col.regions=spectrum, colorkey=TRUE))


###################################################
### code chunk number 14: use
###################################################
use <- 75:274


###################################################
### code chunk number 15: image3 (eval = FALSE)
###################################################
## image(ld.ceph$D.prime[use,use], lwd=0)


###################################################
### code chunk number 16: image3a
###################################################
print(image(ld.ceph$D.prime[use,use], lwd=0))


###################################################
### code chunk number 17: image4 (eval = FALSE)
###################################################
## image(ld.ceph$R.squared[use,use], lwd=0)


###################################################
### code chunk number 18: image4a
###################################################
print(image(ld.ceph$R.squared[use,use], lwd=0))


###################################################
### code chunk number 19: distance
###################################################
pos <- support.ld$Position
diags <- vector("list", 100)
for (i in 1:100) diags[[i]] <- pos[(i+1):603] - pos[1:(603-i)]
dist <- bandSparse(603, k=1:100, diagonals=diags)


###################################################
### code chunk number 20: values
###################################################
distance <- dist@x
D.prime <- ld.ceph$D.prime@x
R.squared <- ld.ceph$R.squared@x


###################################################
### code chunk number 21: drplot
###################################################
plot(hexbin(D.prime^2, R.squared))


###################################################
### code chunk number 22: dpplot1
###################################################
plot(hexbin(distance, D.prime, xbin=10))


###################################################
### code chunk number 23: dpplot2
###################################################
plot(hexbin(distance, R.squared, xbin=10))


###################################################
### code chunk number 24: two
###################################################
snp1 <- as(ceph.1mb[,1], "character")
snp5 <- as(ceph.1mb[,5], "character")
tab33 <- table(snp1, snp5)
tab33


###################################################
### code chunk number 25: twoDR
###################################################
ld.ceph$D.prime[1,5]
ld.ceph$R.squared[1,5]


###################################################
### code chunk number 26: digits
###################################################
options(digits=4)


###################################################
### code chunk number 27: OR
###################################################
OR <- ld(ceph.1mb[,1], ceph.1mb[,5], stats="OR")
OR
AABB <- tab33[2,2]*OR/(1+OR)
ABBA <- tab33[2,2]*1/(1+OR)
AABB
ABBA


###################################################
### code chunk number 28: twoxtwo
###################################################
a <- format(2*tab33[1,1] + tab33[1,2] + tab33[2,1] + AABB, digits=4)
b <- format(2*tab33[1,3] + tab33[1,2] + tab33[2,3] + ABBA, digits=4)
c <- format(2*tab33[3,1] + tab33[2,1] + tab33[3,2] + ABBA, digits=4)
d <- format(2*tab33[3,3] + tab33[3,2] + tab33[2,3] + AABB, digits=4)


###################################################
### code chunk number 29: extent
###################################################
lr100 <- c(68:167, 169:268)
D.prime <- ld(ceph.1mb[,168], ceph.1mb[,lr100], stats="D.prime")
where <- pos[lr100]


###################################################
### code chunk number 30: eplot
###################################################
plot(where, D.prime)
lines(where, smooth(D.prime))


###################################################
### code chunk number 31: ldregion
###################################################
use <- pos>=1.579e7 & pos<=1.587e7
r2 <- forceSymmetric(ld.ceph$R.squared[use, use])


###################################################
### code chunk number 32: dist
###################################################
D <- as.dist(1-r2)
hc <- hclust(D, method="complete")


###################################################
### code chunk number 33: clplot
###################################################
par(cex=0.5)
plot(hc)


###################################################
### code chunk number 34: clusters
###################################################
clusters <- cutree(hc, h=0.5)
head(clusters)
table(clusters)


