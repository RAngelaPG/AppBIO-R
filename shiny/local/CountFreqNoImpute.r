
CounToFreqNI<-function(data_set){
## arranging data
#library(plyr)
#library(reshape2)
#library(reshape)
names(data_set)[1:2] <- c('Allele','SNP')
AlleleID_wohifen <- gsub(pattern = '--', replacement = '-', x = data_set$Allele)
AlleleID_markers <- gsub("([ACTG]>[ACTG]).*", "\\1", AlleleID_wohifen)
## replace collumns
data_set$Allele <- gsub(pattern = '-', replacement = '--', x = AlleleID_markers)
data_set$SNP <- ifelse(is.na(data_set$SNP), 1, 2)
markerID <- 'Allele'
alleleID <- 'SNP'
colsID <- 1:2
print("Prepare long format...")
## reshapes data to long format
long_shape <- melt(data = data_set, id.vars = names(data_set)[colsID],variable.name = 'Gid', value.name = 'Count')  
names(long_shape) <- c('Allele','SNP', 'Gid', 'Count')    
print("Cast alelles...")
## casts w.r.t. alleles
cast_shape <- dcast(long_shape,formula = list(c(markerID,'Gid'), c(alleleID)),var.value='Count')
names(cast_shape) <- c(markerID, 'Gid', 'a1', 'a2')
rm(long_shape)
gc()
print("Calculate allele frequencies...")
## obtains the allele frequencies **NaN as 0**
cast_freq <- mutate(cast_shape,
                    f1 = ifelse(is.nan(a1 / (a1 + a2)), NA, a1 / (a1 + a2)),
                    f2 = ifelse(is.nan(a2 / (a1 + a2)), NA, a2 / (a1 + a2)),
                    a1 = NULL, a2 = NULL) # drop count data
rm(cast_shape)
gc()
ND=cast(cast_freq[,-4],Allele~Gid)
rm(cast_freq)
gc()
names(ND)[1]="AlleleID"
write.csv(ND,"DataforBIO.csv",row.names = F, quote=F)
shinyalert(title = "Important message", 
           text="Successful recodification to allele frequency. \n  You can find the recodification file like 'DataforBIO.csv' ",closeOnEsc = FALSE, 
           type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
)
return(ND)
}