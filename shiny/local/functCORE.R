funchunt=function(datosgen,dir_fileGen,dir_filePhen,dir_fileDist,size1,ENMR,ENCE,ENGD,ENPD,ANMR,ANCE,ANGD,ANPD,EEMR,EECE,EEGD,EEPD,SH,HE,CV){

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

Gen=core1$sel
core1[[1]]=NULL
sumy=unlist(core1)
nmcol=datosgen[,1]
datosgen=as.data.frame(t(datosgen[,-1]))
names(datosgen)=nmcol

if (dir_fileGen[1]!="none" & dir_filePhen[1]=="none" & dir_fileDist[1]=="none"){
  Gen1=datosgen[,as.numeric(Gen)]    
  write.csv(Gen1,paste("CoreSubset_Gen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
  Gen1=datosgen[,as.numeric(Gen)]
  Phen1=pheno.file[as.numeric(Gen),]
    
  write.csv(Gen1,paste("CoreSubset_Gen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(Phen1,paste("CoreSubset_Phen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

if (dir_fileGen[1]!="none" & dir_filePhen[1]!="none" & dir_fileDist[1]!="none"){
  Gen1=datosgen[,as.numeric(Gen)] 
  Phen1=pheno.file[as.numeric(Gen),]
  colnames(dist.file)=c("names",dist.file[,1])
  rownames(dist.file)=dist.file[,1]
  Dist1=dist.file[as.numeric(Gen),as.numeric(Gen)+1]  
  
  write.csv(Gen1,paste("CoreSubset_Gen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(Phen1,paste("CoreSubset_Phen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(Dist1,paste("CoreSubset_Dist.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
  
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]!="none" & dir_fileDist[1]=="none"){
  
  Phen1=pheno.file[as.numeric(Gen),]
    
  write.csv(Phen1,paste("CoreSubset_Phen.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

if (dir_fileGen[1]=="none" & dir_filePhen[1]=="none" & dir_fileDist[1]!="none"){
  colnames(dist.file)=c("names",dist.file[,1])
  rownames(dist.file)=dist.file[,1]
  Dist1=dist.file[as.numeric(Gen),as.numeric(Gen)+1]
 
  write.csv(Dist1,paste("CoreSubset_Dist.csv",sep=""),row.names=TRUE,quote=FALSE)
  write.csv(sumy,paste("SummaryCoreSubset.csv",sep=""),quote=FALSE)
}

}
