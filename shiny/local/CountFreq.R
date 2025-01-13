
CounToFreq<-function(data_set){
#setwd('D:\\MyDocuments\\CIMMYT\\2020\\Petrolli\\PruebaCount')
## arranging data
#data_set <- read.table(file = 'UseData.csv',header = TRUE, sep = ',', dec = '.', as.is = TRUE,na.strings = c('',".",NA))
#library(plyr)
#library(reshape2)
names(data_set)[1:2] <- c('Allele','SNP')
AlleleID_wohifen <- gsub(pattern = '--', replacement = '-', x = data_set$Allele)
AlleleID_markers <- gsub("([ACTG]>[ACTG]).*", "\\1", AlleleID_wohifen)
## replace collumns
data_set$Allele <- gsub(pattern = '-', replacement = '--', x = AlleleID_markers)
data_set$SNP <- ifelse(is.na(data_set$SNP), 1, 2)
markerID <- 'Allele'
alleleID <- 'SNP'
colsID <- 1:2
maf <- .01
## reshapes data to long format
long_shape <- melt(data = data_set, id.vars = names(data_set)[colsID],variable.name = 'Gid', value.name = 'Count')  
## counts alleles & markers
counts_allele <- ddply(.data = long_shape, .variables = c(markerID, alleleID), .fun = summarize, by_allele = sum(Count))
counts_marker <- ddply(.data = counts_allele, .variables = c(markerID), .fun = summarize, by_marker = sum(by_allele))
counts_joint <- merge(counts_allele, counts_marker, by = markerID)
## casts w.r.t. alleles
cast_shape <- dcast(long_shape,formula = list(c(markerID,'Gid'), c(alleleID)),value.var = 'Count')
names(cast_shape) <- c(markerID, 'Gid', 'a1', 'a2')    
## obtains coverage statistics
counts_cover <- ddply(.data = cast_shape, .variables = c(markerID),.fun = summarize, coverage = sum((a1 + a2) > 0))
counts_all <- merge(counts_joint, counts_cover, by = markerID)
counts_stats <- mutate(counts_all,aveMaf = by_allele / by_marker,meanCover = by_marker / (2 * coverage),aveCover = aveMaf * meanCover)    
## filters w.r.t. MAF
maf_reached <- subset(counts_stats, ( aveMaf > maf ) & ( aveMaf < (1 - maf) ))
subset_maf <- long_shape[long_shape[[markerID]] %in% unique(maf_reached[[markerID]]),]
## filters w.r.t. missings
missing_all <- melt(subset(cast_shape, (a1 + a2) == 0 ),id.vars = c(markerID, 'Gid'),variable.name = alleleID, value.name = 'void')
levels(missing_all[[alleleID]]) <- 1:2       
## reduces statistics to missings
missing_stats <- merge(missing_all, counts_stats, by = c(markerID, alleleID))    
## catches markers to be imputed
missing_impute <- merge(subset_maf, missing_stats,by = c(markerID, alleleID, 'Gid'), all.x = TRUE)
## imputation
missing_impute$imputed <- with(missing_impute,ifelse(test = is.na(aveCover), yes = Count, no = aveCover))    
## casts to original shape
imputed <- dcast(missing_impute, list(c(markerID, alleleID), c('Gid')),value.var = 'imputed')
long_shape <- melt(imputed, id.vars = c('SNP', 'Allele'),variable.name = 'Gid', value.name = 'Count')
names(long_shape)=c("Allele","SNP","Gid","Count")
## casts w.r.t Alleles
cast_shape <- dcast(long_shape, formula = SNP + Gid ~ Allele,value.var = 'Count')
names(cast_shape) <- c('SNP', 'Gid', 'a1', 'a2')
## obtains the allele frequencies **NaN as 0**
cast_freq <- mutate(cast_shape,
                    f1 = ifelse(is.nan(a1 / (a1 + a2)), NA, a1 / (a1 + a2)),
                    f2 = ifelse(is.nan(a2 / (a1 + a2)), NA, a2 / (a1 + a2)),
                    a1 = NULL, a2 = NULL) # drop count data
## reshapes the casted frequencies
long_freq <- melt(cast_freq, id.vars = c('SNP', 'Gid'),variable.name = 'Allele', value.name = 'freq')
levels(long_freq$Allele) <- 1:2
long_freq$Allele <- as.numeric(levels(long_freq$Allele)[long_freq$Allele])
## casts w.r.t "gids" -- reverses to original shape
freqs <- dcast(long_freq, formula = SNP + Allele ~ Gid,value.var = 'freq')
freqs <- arrange(freqs, SNP, Allele) # arrange output
ND=freqs[which(freqs[,2]==2),-2]
names(ND)[1]="AlleleID"
write.csv(ND,"DataforBIO.csv",row.names = F, quote=F)
shinyalert(title = "Important message", 
           text="Successful recodification to allele frequency. \n  You can find the recodification file like 'DataforBIO.csv' ",closeOnEsc = FALSE, 
           type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
)
return(ND)
}