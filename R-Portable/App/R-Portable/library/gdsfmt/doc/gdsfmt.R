## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("gdsfmt")

## ----eval=FALSE---------------------------------------------------------------
#  library("devtools")
#  install_github("zhengxwen/gdsfmt")

## -----------------------------------------------------------------------------
library(gdsfmt)
gfile <- createfn.gds("test.gds")

## ----echo=FALSE---------------------------------------------------------------
set.seed(1000)

## -----------------------------------------------------------------------------
add.gdsn(gfile, "int", val=1:10000)
add.gdsn(gfile, "double", val=seq(1, 1000, 0.4))
add.gdsn(gfile, "character", val=c("int", "double", "logical", "factor"))
add.gdsn(gfile, "logical", val=rep(c(TRUE, FALSE, NA), 50), visible=FALSE)
add.gdsn(gfile, "factor", val=as.factor(c(NA, "AA", "CC")), visible=FALSE)
add.gdsn(gfile, "bit2", val=sample(0:3, 1000, replace=TRUE), storage="bit2")

# list and data.frame
add.gdsn(gfile, "list", val=list(X=1:10, Y=seq(1, 10, 0.25)))
add.gdsn(gfile, "data.frame", val=data.frame(X=1:19, Y=seq(1, 10, 0.5)))

## -----------------------------------------------------------------------------
folder <- addfolder.gdsn(gfile, "folder")
add.gdsn(folder, "int", val=1:1000)
add.gdsn(folder, "double", val=seq(1, 100, 0.4), visible=FALSE)

## -----------------------------------------------------------------------------
gfile

## -----------------------------------------------------------------------------
print(gfile, all=TRUE)

## -----------------------------------------------------------------------------
index.gdsn(gfile, "int")
index.gdsn(gfile, "list/Y")
index.gdsn(gfile, "folder/int")

## -----------------------------------------------------------------------------
# close the GDS file
closefn.gds(gfile)

## -----------------------------------------------------------------------------
gfile <- createfn.gds("test.gds")

## -----------------------------------------------------------------------------
n <- add.gdsn(gfile, "I1", val=matrix(1:15, nrow=3))
show(n)

## -----------------------------------------------------------------------------
write.gdsn(n, rep(0,5), start=c(2,1), count=c(1,-1))
show(n)

## -----------------------------------------------------------------------------
append.gdsn(n, 16:24)
show(n)

## -----------------------------------------------------------------------------
# initialize
n <- add.gdsn(gfile, "mat", matrix(1:48, 6))
show(n)

# substitute
assign.gdsn(n, .value=c(9:14,35:40), .substitute=NA)
show(n)

# subset
assign.gdsn(n, seldim=list(rep(c(TRUE, FALSE),3), rep(c(FALSE, TRUE),4)))
show(n)

# initialize and subset
n <- add.gdsn(gfile, "mat", matrix(1:48, 6), replace=TRUE)
assign.gdsn(n, seldim=list(c(4,2,6,NA), c(5,6,NA,2,8,NA,4)))
show(n)

# initialize and sort into descending order
n <- add.gdsn(gfile, "mat", matrix(1:48, 6), replace=TRUE)
assign.gdsn(n, seldim=list(6:1, 8:1))
show(n)

## -----------------------------------------------------------------------------
(n2 <- add.gdsn(gfile, "I2", storage="int", valdim=c(100, 2000)))

for (i in 1:2000)
{
    write.gdsn(n2, seq.int(100*(i-1)+1, length.out=100),
        start=c(1,i), count=c(-1,1))
}

show(n2)

## -----------------------------------------------------------------------------
(n3 <- add.gdsn(gfile, "I3", storage="int", valdim=c(100, 0), compress="ZIP"))

for (i in 1:2000)
{
    append.gdsn(n3, seq.int(100*(i-1)+1, length.out=100))
}

readmode.gdsn(n3)  # finish writing with the compression algorithm
show(n3)

## -----------------------------------------------------------------------------
# close the GDS file
closefn.gds(gfile)

## -----------------------------------------------------------------------------
gfile <- createfn.gds("test.gds")
add.gdsn(gfile, "I1", val=matrix(1:20, nrow=4))
add.gdsn(gfile, "I2", val=1:100)
closefn.gds(gfile)

## -----------------------------------------------------------------------------
gfile <- openfn.gds("test.gds")
n <- index.gdsn(gfile, "I1")

read.gdsn(n)

## -----------------------------------------------------------------------------
# read a subset
read.gdsn(n, start=c(2, 2), count=c(2, 3))

read.gdsn(n, start=c(2, 2), count=c(2, 3), .value=c(6,15), .substitute=NA)

## -----------------------------------------------------------------------------
# read a subset
readex.gdsn(n, list(c(FALSE,TRUE,TRUE,FALSE), c(TRUE,FALSE,TRUE,FALSE,TRUE)))

readex.gdsn(n, list(c(1,4,3,NA), c(2,NA,3,1,3,1)))

readex.gdsn(n, list(c(1,4,3,NA), c(2,NA,3,1,3,1)), .value=NA, .substitute=-1)

## -----------------------------------------------------------------------------
apply.gdsn(n, margin=1, FUN=print, as.is="none")

apply.gdsn(n, margin=2, FUN=print, as.is="none")

# close the GDS file
closefn.gds(gfile)

## -----------------------------------------------------------------------------
gfile <- createfn.gds("test.gds")
n1 <- add.gdsn(gfile, "I1", val=1:100)
n2 <- add.gdsn(gfile, "I2", val=matrix(1:20, nrow=4))

gfile

## -----------------------------------------------------------------------------
fout <- file("text.txt", "wt")
apply.gdsn(n1, 1, FUN=cat, as.is="none", file=fout, fill=TRUE)
close(fout)

scan("text.txt")

## -----------------------------------------------------------------------------
fout <- file("text.txt", "wt")
apply.gdsn(n2, 1, FUN=cat, as.is="none", file=fout, fill=4194304)
close(fout)

readLines("text.txt")

## -----------------------------------------------------------------------------
n.t <- add.gdsn(gfile, "transpose", storage="int", valdim=c(5,0))

# apply the function over rows of matrix
apply.gdsn(n2, margin=1, FUN=c, as.is="gdsnode", target.node=n.t)

# matrix transpose
read.gdsn(n.t)

# close the GDS file
closefn.gds(gfile)

## -----------------------------------------------------------------------------
set.seed(1000)
val <- sample(seq(0,1,0.001), 50000, replace=TRUE)
head(val)

gfile <- createfn.gds("test.gds")

add.gdsn(gfile, "N1", val=val)
add.gdsn(gfile, "N2", val=val, compress="ZIP", closezip=TRUE)
add.gdsn(gfile, "N3", val=val, storage="float")
add.gdsn(gfile, "N4", val=val, storage="float", compress="ZIP", closezip=TRUE)
add.gdsn(gfile, "N5", val=val, storage="packedreal16", scale=0.001, offset=0)
add.gdsn(gfile, "N6", val=val, storage="packedreal16", scale=0.001, offset=0, compress="ZIP", closezip=TRUE)

gfile

## ----echo=FALSE---------------------------------------------------------------
KB <- function(i)
{
    s <- objdesp.gdsn(index.gdsn(gfile, paste("N", i, sep="")))$size
    sprintf("%.1f KB", s/1000)
}

Ratio <- function(i)
{
    s <- objdesp.gdsn(index.gdsn(gfile, paste("N", i, sep="")))$size
    r <- 100*s / (8*length(val))
    sprintf("%.1f%%", r)
}

Epsilon <- function(i)
{
    ans <- mean(abs(val - read.gdsn(index.gdsn(gfile, paste0("N",i)))))
    sprintf("%0.3g", ans)
}

## -----------------------------------------------------------------------------
# close the GDS file
closefn.gds(gfile)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(100)
#  # 10,000,000 random 0,1 sequence of 32-bit integers
#  val <- sample.int(2, 10*1000*1000, replace=TRUE) - 1L
#  table(val)

## ----eval=FALSE---------------------------------------------------------------
#  # cteate a GDS file
#  f <- createfn.gds("test.gds")
#  
#  # compression algorithms (LZMA_ra:32K is the lower bound of LZMA_ra)
#  compression <- c("", "ZIP.max", "ZIP_ra.max:16K", "LZ4.max", "LZ4_ra.max:16K", "LZMA", "LZMA_ra:32K")
#  
#  # save
#  for (i in 1:length(compression))
#  	print(add.gdsn(f, paste0("I", i), val=val, compress=compression[i], closezip=TRUE))
#  
#  # close the file
#  closefn.gds(f)
#  
#  cleanup.gds("test.gds")

## ----eval=FALSE---------------------------------------------------------------
#  # open the GDS file
#  f <- openfn.gds("test.gds")
#  
#  # 10,000 random positions
#  set.seed(1000)
#  idx <- sample.int(length(val), 10000)
#  
#  # enumerate each compression method
#  dat <- vector("list", length(compression))
#  for (i in seq_len(length(compression)))
#  {
#  	cat("Compression:", compression[i], "\n")
#  	n <- index.gdsn(f, paste0("I", i))
#  	print(system.time({
#  		dat[[i]] <- sapply(idx, FUN=function(k) read.gdsn(n, start=k, count=1L))
#  	}))
#  }
#  
#  # check
#  for (i in seq_len(length(compression)))
#  	stopifnot(identical(dat[[i]], dat[[1L]]))
#  
#  # close the file
#  closefn.gds(f)

## -----------------------------------------------------------------------------
# create a GDS file
f <- createfn.gds("test.gds")
set.seed(1000)
m <- matrix(sample(c(0:2), 56, replace=T), nrow=4)
(n <- add.gdsn(f, "sparse", m, storage="sp.int"))
# get a dgCMatrix sparse matrix (.sparse=TRUE by default)
read.gdsn(n)
# get a dense matrix
read.gdsn(n, .sparse=FALSE)
closefn.gds(f)

## -----------------------------------------------------------------------------
# create a GDS file
f <- createfn.gds("test.gds")
n <- add.gdsn(f, "raw", rnorm(1115), compress="ZIP", closezip=TRUE)
digest.gdsn(n, action="add")
print(f, attribute=TRUE)
closefn.gds(f)

## -----------------------------------------------------------------------------
f <- openfn.gds("test.gds")
n <- index.gdsn(f, "raw")

get.attr.gdsn(n)$md5
digest.gdsn(n, action="verify")  # NA indicates "not applicable"

closefn.gds(f)

## -----------------------------------------------------------------------------
sessionInfo()

## ----echo=FALSE---------------------------------------------------------------
unlink(c("test.gds", "text.txt"), force=TRUE)

