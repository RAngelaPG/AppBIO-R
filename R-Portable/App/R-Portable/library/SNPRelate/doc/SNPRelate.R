## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("gdsfmt")
#  BiocManager::install("SNPRelate")

## ----eval=FALSE---------------------------------------------------------------
#  library("devtools")
#  install_github("zhengxwen/gdsfmt")
#  install_github("zhengxwen/SNPRelate")

## -----------------------------------------------------------------------------
# Load the R packages: gdsfmt and SNPRelate
library(gdsfmt)
library(SNPRelate)

## -----------------------------------------------------------------------------
snpgdsSummary(snpgdsExampleFileName())

## -----------------------------------------------------------------------------
# Open a GDS file
(genofile <- snpgdsOpen(snpgdsExampleFileName()))

## -----------------------------------------------------------------------------
# Get the attributes of chromosome coding
get.attr.gdsn(index.gdsn(genofile, "snp.chromosome"))

## -----------------------------------------------------------------------------
# Take out genotype data for the first 3 samples and the first 5 SNPs
(g <- read.gdsn(index.gdsn(genofile, "genotype"), start=c(1,1), count=c(5,3)))

## ----eval=FALSE---------------------------------------------------------------
#  g <- snpgdsGetGeno(genofile, sample.id=..., snp.id=...)

## -----------------------------------------------------------------------------
# Get the attribute of genotype
get.attr.gdsn(index.gdsn(genofile, "genotype"))

## -----------------------------------------------------------------------------
# Take out snp.id
head(read.gdsn(index.gdsn(genofile, "snp.id")))
# Take out snp.rs.id
head(read.gdsn(index.gdsn(genofile, "snp.rs.id")))

## -----------------------------------------------------------------------------
# Read population information
pop <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))
table(pop)

# Close the GDS file
snpgdsClose(genofile)

## -----------------------------------------------------------------------------
# Load data
data(hapmap_geno)

# Create a gds file
snpgdsCreateGeno("test.gds", genmat = hapmap_geno$genotype,
    sample.id = hapmap_geno$sample.id, snp.id = hapmap_geno$snp.id,
    snp.chromosome = hapmap_geno$snp.chromosome,
    snp.position = hapmap_geno$snp.position,
    snp.allele = hapmap_geno$snp.allele, snpfirstdim=TRUE)

# Open the GDS file
(genofile <- snpgdsOpen("test.gds"))

# Close the GDS file
snpgdsClose(genofile)

## ----eval=FALSE---------------------------------------------------------------
#  # Create a new GDS file
#  newfile <- createfn.gds("your_gds_file.gds")
#  
#  # add a flag
#  put.attr.gdsn(newfile$root, "FileFormat", "SNP_ARRAY")
#  
#  # Add variables
#  add.gdsn(newfile, "sample.id", sample.id)
#  add.gdsn(newfile, "snp.id", snp.id)
#  add.gdsn(newfile, "snp.chromosome", snp.chromosome)
#  add.gdsn(newfile, "snp.position", snp.position)
#  add.gdsn(newfile, "snp.allele", c("A/G", "T/C", ...))
#  
#  #####################################################################
#  # Create a snp-by-sample genotype matrix
#  
#  # Add genotypes
#  var.geno <- add.gdsn(newfile, "genotype",
#      valdim=c(length(snp.id), length(sample.id)), storage="bit2")
#  
#  # Indicate the SNP matrix is snp-by-sample
#  put.attr.gdsn(var.geno, "snp.order")
#  
#  # Write SNPs into the file sample by sample
#  for (i in 1:length(sample.id))
#  {
#      g <- ...
#      write.gdsn(var.geno, g, start=c(1,i), count=c(-1,1))
#  }
#  
#  #####################################################################
#  # OR, create a sample-by-snp genotype matrix
#  
#  # Add genotypes
#  var.geno <- add.gdsn(newfile, "genotype",
#      valdim=c(length(sample.id), length(snp.id)), storage="bit2")
#  
#  # Indicate the SNP matrix is sample-by-snp
#  put.attr.gdsn(var.geno, "sample.order")
#  
#  # Write SNPs into the file sample by sample
#  for (i in 1:length(snp.id))
#  {
#      g <- ...
#      write.gdsn(var.geno, g, start=c(1,i), count=c(-1,1))
#  }
#  
#  # Get a description of chromosome codes
#  #   allowing to define a new chromosome code, e.g., snpgdsOption(Z=27)
#  option <- snpgdsOption()
#  var.chr <- index.gdsn(newfile, "snp.chromosome")
#  put.attr.gdsn(var.chr, "autosome.start", option$autosome.start)
#  put.attr.gdsn(var.chr, "autosome.end", option$autosome.end)
#  for (i in 1:length(option$chromosome.code))
#  {
#      put.attr.gdsn(var.chr, names(option$chromosome.code)[i],
#          option$chromosome.code[[i]])
#  }
#  
#  # Add your sample annotation
#  samp.annot <- data.frame(sex = c("male", "male", "female", ...),
#      pop.group = c("CEU", "CEU", "JPT", ...), ...)
#  add.gdsn(newfile, "sample.annot", samp.annot)
#  
#  # Add your SNP annotation
#  snp.annot <- data.frame(pass=c(TRUE, TRUE, FALSE, FALSE, TRUE, ...), ...)
#  add.gdsn(newfile, "snp.annot", snp.annot)
#  
#  # Close the GDS file
#  closefn.gds(newfile)

## -----------------------------------------------------------------------------
# The PLINK BED file, using the example in the SNPRelate package
bed.fn <- system.file("extdata", "plinkhapmap.bed.gz", package="SNPRelate")
fam.fn <- system.file("extdata", "plinkhapmap.fam.gz", package="SNPRelate")
bim.fn <- system.file("extdata", "plinkhapmap.bim.gz", package="SNPRelate")

## ----eval=FALSE---------------------------------------------------------------
#  bed.fn <- "C:/your_folder/your_plink_file.bed"
#  fam.fn <- "C:/your_folder/your_plink_file.fam"
#  bim.fn <- "C:/your_folder/your_plink_file.bim"

## -----------------------------------------------------------------------------
# Convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "test.gds")

# Summary
snpgdsSummary("test.gds")

## -----------------------------------------------------------------------------
# The VCF file, using the example in the SNPRelate package
vcf.fn <- system.file("extdata", "sequence.vcf", package="SNPRelate")

## ----eval=FALSE---------------------------------------------------------------
#  vcf.fn <- "C:/your_folder/your_vcf_file.vcf"

## -----------------------------------------------------------------------------
# Reformat
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

# Summary
snpgdsSummary("test.gds")

## -----------------------------------------------------------------------------
# Open the GDS file
genofile <- snpgdsOpen(snpgdsExampleFileName())

## -----------------------------------------------------------------------------
# Get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))

table(pop_code)

# Display the first six values
head(pop_code)

## -----------------------------------------------------------------------------
set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
str(snpset)
names(snpset)

# Get all selected snp id
snpset.id <- unlist(unname(snpset))
head(snpset.id)

## -----------------------------------------------------------------------------
# Run PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

## -----------------------------------------------------------------------------
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

# Draw
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

## -----------------------------------------------------------------------------
# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# assume the order of sample IDs is as the same as population codes
head(cbind(sample.id, pop_code))

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# Make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

# Draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)

## ----fig.width=8, fig.height=5, fig.align='center'----------------------------
library(MASS)

datpop <- factor(pop_code)[match(pca$sample.id, sample.id)]
parcoord(pca$eigenvect[,1:16], col=datpop)

## ----fig.width=8, fig.height=4, fig.align='center'----------------------------
# Get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)

savepar <- par(mfrow=c(2,1), mai=c(0.45, 0.55, 0.1, 0.25))
for (i in 1:2)
{
    plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i),
        col=chr, pch="+")
}
par(savepar)

## -----------------------------------------------------------------------------
# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# Two populations: HCB and JPT
flag <- pop_code %in% c("HCB", "JPT")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
v <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel),
    method="W&C84")
v$Fst        # Weir and Cockerham weighted Fst estimate
v$MeanFst    # Weir and Cockerham mean Fst estimate
summary(v$FstSNP)

# Multiple populations: CEU HCB JPT YRI
#   we should remove offsprings
father <- read.gdsn(index.gdsn(genofile, "sample.annot/father.id"))
mother <- read.gdsn(index.gdsn(genofile, "sample.annot/mother.id"))
flag <- (father=="") & (mother=="")
samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
v <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel),
    method="W&C84")
v$Fst        # Weir and Cockerham weighted Fst estimate
v$MeanFst    # Weir and Cockerham mean Fst estimate
summary(v$FstSNP)

## -----------------------------------------------------------------------------
# YRI samples
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
YRI.id <- sample.id[pop_code == "YRI"]

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=YRI.id, snp.id=snpset.id,
    maf=0.05, missing.rate=0.05, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
    xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# Estimate IBD coefficients
set.seed(100)
snp.id <- sample(snpset.id, 1500)  # random 1500 SNPs
ibd <- snpgdsIBDMLE(genofile, sample.id=YRI.id, snp.id=snp.id,
    maf=0.05, missing.rate=0.05, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
    xlab="k0", ylab="k1", main="YRI samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# Incorporate with pedigree information
family.id <- read.gdsn(index.gdsn(genofile, "sample.annot/family.id"))
family.id <- family.id[match(YRI.id, sample.id)]
table(family.id)

ibd.robust <- snpgdsIBDKING(genofile, sample.id=YRI.id,
    family.id=family.id, num.thread=2)
names(ibd.robust)

# Pairs of individuals
dat <- snpgdsIBDSelection(ibd.robust)
head(dat)

plot(dat$IBS0, dat$kinship, xlab="Proportion of Zero IBS",
    ylab="Estimated Kinship Coefficient (KING-robust)")

## -----------------------------------------------------------------------------
ibs <- snpgdsIBS(genofile, num.thread=2)

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# individulas in the same population are clustered together
pop.idx <- order(pop_code)

image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code)

plot(x, y, col=race, xlab = "", ylab = "",
    main = "Multidimensional Scaling Analysis (IBS)")
legend("topleft", legend=levels(race), pch="o", text.col=1:nlevels(race))

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2))

# Determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none", main="HapMap Phase II")
table(rv$samp.group)

## ----fig.width=5, fig.height=5, fig.align='center'----------------------------
# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

plot(rv2$dendrogram, leaflab="none", main="HapMap Phase II")
legend("topright", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)

# Close the GDS file
snpgdsClose(genofile)

## -----------------------------------------------------------------------------
sessionInfo()

## ----echo=FALSE---------------------------------------------------------------
unlink("test.gds", force=TRUE)

