report_stats<-function(objectS){
	pA <- colMeans(objectS, na.rm = TRUE) / 2 
	He_loci <- 2 * pA * (1 - pA)
	mean_He <- round(mean(He_loci, na.rm = TRUE),3)
	sd_He <- round(sd(He_loci, na.rm = TRUE),3)
	sd_val <- sd(He_loci, na.rm = TRUE)
	mean_val <- mean(He_loci, na.rm = TRUE)
	CV<-round((sd_val / mean_val) * 100,3)
	het_obs <- apply(objectS, 1, function(ind) mean(ind == 1, na.rm = TRUE))  # 1 = heterocigoto
	mean_Ho<-round(mean(het_obs, na.rm= TRUE),3)
	sd_Ho<-round(sd(het_obs, na.rm= TRUE),3)
	shannon_index <- function(p) {
	p <- p[p > 0]  # eliminar ceros para evitar log(0)
	-sum(p * log(p))
	}
	shannon_values <- numeric(length(pA))
	for (i in seq_along(pA)) {
	freqs <- c(pA[i], 1 - pA[i])  # alelo alt y ref
	shannon_values[i] <- shannon_index(freqs)
	}
	mean_shannon <- round(mean(shannon_values, na.rm = TRUE),3)
	sd_shannon <- round(sd(shannon_values, na.rm = TRUE),3)
	return(c(mean_He,sd_He,mean_Ho,sd_Ho,CV,mean_shannon,sd_shannon))
}

funchunt<-function(datosgen,dir_fileGen,dir_filePhen,dir_fileDist,size1,ENMR,ENCE,ENGD,ENPD,ANMR,ANCE,ANGD,ANPD,EEMR,EECE,EEGD,EEPD,SH,HE,CV){
#save(datosgen,dir_fileGen,dir_filePhen,dir_fileDist,size1,ENMR,ENCE,ENGD,ENPD,ANMR,ANCE,ANGD,ANPD,EEMR,EECE,EEGD,EEPD,SH,HE,CV,file="testCH.RData")

if (dir_fileGen[1]!="none" & dir_filePhen[1]=="none" & dir_fileDist[1]=="none"){
   geno.file=datosgen
   my.data <- genotypes(geno.file,format="frequency")
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
    geno.file=datosgen
    pheno.file=read.autodelim(dir_filePhen$datapath)
    my.data <- coreHunterData(genotypes(geno.file,format="frequency"),phenotypes(pheno.file))
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]!="none"){
   geno.file=datosgen
   pheno.file=read.autodelim(dir_filePhen$datapath)
   dist.file=read.autodelim(dir_fileDist$datapath)
   my.data <- coreHunterData(genotypes(geno.file,format="frequency"),phenotypes(pheno.file),distances(dist.file))
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
   pheno.file=read.autodelim(dir_filePhen$datapath)
   my.data <- phenotypes(pheno.file)
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]=="none" & dir_fileDist[1]!="none"){
   dist.file=read.autodelim(dir_fileDist$datapath)
   my.data <- distances(dist.file)
}

print("Optimizing...")
values=as.numeric(c(ENMR,ENCE,ENGD,ENPD,ANMR,ANCE,ANGD,ANPD,EEMR,EECE,EEGD,EEPD,SH,HE,CV))
objall=list(objective("EN", "MR", weight=as.numeric(ENMR)),
            objective("EN", "CE", weight=as.numeric(ENCE)),
            objective("EN", "GD", weight=as.numeric(ENGD)),
            objective("EN", "PD", weight=as.numeric(ENPD)),
            objective("AN", "MR", weight=as.numeric(ANMR)),
            objective("AN", "CE", weight=as.numeric(ANCE)),
            objective("AN", "GD", weight=as.numeric(ANGD)),
            objective("AN", "PD", weight=as.numeric(ANPD)),
            objective("EE", "MR", weight=as.numeric(EEMR)),
            objective("EE", "CE", weight=as.numeric(EECE)),
            objective("EE", "GD", weight=as.numeric(EEGD)),
            objective("EE", "PD", weight=as.numeric(EEPD)),
            objective("SH", weight=as.numeric(SH)),
            objective("HE", weight=as.numeric(HE)),
            objective("CV", weight=as.numeric(CV)))
objuse=list()
use=which(values!=0)
for(i in 1:length(use)){
objuse[[i]]=objall[[use[i]]]
}

core1=sampleCore(my.data,obj=objuse,size=as.numeric(size1),indices=TRUE)

matdist<-function(datos,metric){
	if (metric=="Rogers"){
		nacc=dim(datos)[2]
		fr=as.matrix(datos)												## recover the marker information
		frn=fr
		frn[!is.na(frn)]=1																## no missing values convert to 1
		frn[is.na(frn)]=0																  ## missing values convert to 0
		N=2*crossprod(frn)																## create square matrix markers information
		rm(frn)
		aux=matrix(0,nacc,nacc)
		for(i in 1:(nacc-2)){
			aux[i,]=t(t(c(rep(0,i),sqrt((2*apply((fr[,i]-fr[,-c(1:i)])^2,2,function(x) sum(x,na.rm=T)))/N[i,-c(1:i)]))))
		}
		aux[nacc-1,]=t(t(c(rep(0,nacc-1),sqrt((2*sum((fr[,nacc-1]-fr[,-c(1:nacc-1)])^2,na.rm=T))/N[nacc-1,-c(1:(nacc-1))]))))
		mrdMAT=as.matrix(aux+t(aux))
		colnames(mrdMAT)=colnames(datos)
		rownames(mrdMAT)=colnames(datos)
	}else{		
		rownames(datos)<-datos[,1]
		datos<-datos[,-1]		
		mrdMAT=cluster::daisy(datos,metric="gower")	
		mrdMAT=as.matrix(mrdMAT)		
	}
	return(mrdMAT)
}

Gen=core1$sel
core1[[1]]=NULL
sumy=unlist(core1)
nmcol=datosgen[,1]
datosgen=as.data.frame(t(datosgen[,-1]))
names(datosgen)=nmcol


if (dir_fileGen[1]!="none" & dir_filePhen[1]=="none" & dir_fileDist[1]=="none"){
  print("Calculating MDS...")
  geno_distAll<-matdist(datosgen,"Rogers")  
  Gen1=datosgen[,as.numeric(Gen)]   
  Gen1=as.data.frame(t(t(colnames(Gen1)))) 
  names(Gen1)=c("Gen")    
  #g_full_sub <- geno_distAll[Gen1[,1], Gen1[,1]]    
  #geno_distSub <- matdist(datosgen[,Gen1[,1]],"Rogers")
  #mantel_res <- vegan::mantel(as.dist(g_full_sub), as.dist(geno_distSub), permutations = 999)  
  #mantel_res<-c(mantel_res$statistic,mantel_res$signif)	  
  print("Saving files...")  
  write.csv(Gen1,paste("CoreSubset_Gen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
  print("Calculating MDS...")
  geno_distAll<-matdist(datosgen,"Rogers")
  Gen1=datosgen[,as.numeric(Gen)]
  Gen1=as.data.frame(t(t(colnames(Gen1)))) 
  names(Gen1)=c("Gen")  
  #Phen1=pheno.file[as.numeric(Gen),]
  print("Saving files...")  
  write.csv(Gen1,paste("CoreSubset_GenPhen.csv",sep=""),row.names=TRUE,quote=FALSE)
  #write.csv(Phen1,paste("CoreSubset_Phen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]!="none"){  
  print("Calculating MDS...")
  geno_distAll<-dist.file
  colnames(geno_distAll)=geno_distAll[,1]
  rownames(geno_distAll)=geno_distAll[,1]
  geno_distAll<-as.matrix(geno_distAll[,-1])
  Gen1=datosgen[,as.numeric(Gen)] 
  Gen1=as.data.frame(t(t(colnames(Gen1))))  
  names(Gen1)=c("Gen")
  #Phen1=pheno.file[as.numeric(Gen),]
  #colnames(dist.file)=c("names",dist.file[,1])
  #rownames(dist.file)=dist.file[,1]
  #Dist1=dist.file[as.numeric(Gen),as.numeric(Gen)+1]  
  print("Saving files...")
  write.csv(Gen1,paste("CoreSubset_GenPhenDist.csv",sep=""),row.names=TRUE,quote=FALSE)
  #write.csv(Phen1,paste("CoreSubset_Phen.csv",sep=""),row.names=TRUE,quote=FALSE)
  #write.csv(Dist1,paste("CoreSubset_Dist.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)  
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
  print("Calculating MDS...")
  geno_distAll<-matdist(pheno.file,"gower")
  Phen1=pheno.file[as.numeric(Gen),]
  Gen1=as.data.frame(Phen1[,1])
  names(Gen1)=c("Gen")
   print("Saving files...") 
  write.csv(Gen1,paste("CoreSubset_Phen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]=="none" & dir_fileDist[1]!="none"){
  print("Calculating MDS...")
  geno_distAll<-dist.file
  colnames(geno_distAll)=geno_distAll[,1]
  rownames(geno_distAll)=geno_distAll[,1]
  geno_distAll<-as.matrix(geno_distAll[,-1])
  colnames(dist.file)=c("names",dist.file[,1])
  rownames(dist.file)=dist.file[,1]
  Dist1=dist.file[as.numeric(Gen),as.numeric(Gen)+1]
  Gen1=as.data.frame(rownames(Dist1))
  names(Gen1)=c("Gen")    
  print("Saving files...")
  write.csv(Gen1,paste("CoreSubset_Dist.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

res<-list(Gen1,geno_distAll)
#save(res,file="TDmds.RData")
return(res)
}
