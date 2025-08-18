CounToFreqNI<-function(dt_ff,Allele,newcolnames){
		AlleleID_wohifen <- gsub(pattern = '--', replacement = '-', x = Allele$AlleleID)
		AlleleID_markers <- gsub("([ACTG]>[ACTG]).*", "\\1", AlleleID_wohifen)
		## replace collumns
		Allele <- gsub(pattern = '-', replacement = '--', x = AlleleID_markers)
		rm(AlleleID_wohifen,AlleleID_markers)
		gc()
		SNP=as.numeric(rep(c(NA,2),nrow(dt_ff)/2))
		SNP <- ifelse(is.na(SNP), 1, 2)
		markerID <- 'Allele'
		alleleID <- 'SNP'
		calcF<-function(x){
		long_shape <- data.frame(cbind(Allele,SNP,names(x),as.data.frame(as.ffdf(x))))
		names(long_shape) <- c('Allele','SNP', 'Gid', 'Count')    
		long_shape$Count=as.numeric(long_shape$Count) 
		## casts w.r.t. alleles
		cast_shape <- dcast(long_shape,formula = list(c(markerID,'Gid'), c(alleleID)),var.value='Count')
		names(cast_shape) <- c(markerID, 'Gid', 'a1', 'a2')
		cast_freq <- mutate(cast_shape,
							f1 = ifelse(is.nan(a1 / (a1 + a2)), NA, a1 / (a1 + a2)),
							f2 = ifelse(is.nan(a2 / (a1 + a2)), NA, a2 / (a1 + a2)),
							a1 = NULL, a2 = NULL) # drop count data
		cast_freq=cast_freq[,c(1,3)]
		names(cast_freq)=c("Allele",names(x))
		return(cast_freq[,1:2])
		}
		
		dt_use=ffcolapply(calcF(dt_ff[,i1:i2,drop=FALSE]), X=dt_ff, RETURN=TRUE, CFUN="ccbind", BATCHSIZE=1,VERBOSE=T)
	    dt_use=dt_use[,-c(2,3,seq(3,dim(dt_use)[2],2))]
		rm(dt_ff,SNP)
		gc()
		colnames(dt_use)=newcolnames
		#dt_use=data.frame(cbind(AlleleID=unique(Allele),data.frame(dt_use[,2:dim(dt_use)[2]])))
		fwrite(dt_use,"DataforBIO.csv",row.names = F, quote=F)
		rm(Allele,newcolnames)
		gc()
		shinyalert(title = "Important message", 
				text="Successful recodification to allele frequency. \n  You can find the recodification file like 'DataforBIO.csv' ",closeOnEsc = FALSE, 
				type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
		)
		return(dt_use)

}
