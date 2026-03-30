Biodv=function(file_name,datos,nall,distk,mayorque,menorque,missval,typedata,ht1,ht2,ht3,missvalG,mixture,gap,methodgap,nclust){
#save(file_name,datos,nall,distk,mayorque,menorque,missval,typedata,ht1,ht2,ht3,gap,file="check.RData")

#######################################################
id=as.data.frame(datos[,1])
colnames(id)=c("Markers")
datos=datos[,2:ncol(datos)]
oriMark=nrow(datos)
oriGen=ncol(datos)
#Add for markers and genotypes replicated
bygen=as.data.frame(apply(datos,2,paste,collapse = "" ))

bygen[,1]=as.character(bygen[,1])

idG=as.data.frame(colnames(datos))
colnames(idG)=c("Genotypes")
dupg = duplicated(bygen) | duplicated(bygen, fromLast = T)
dupg=split(which(dupg), bygen[dupg,])
dupg=lapply(dupg,function(x)idG[x,])
if (length(dupg)!=0){
	max.lengthg <- max(sapply(dupg, length))
	dupg <- lapply(dupg, function(v) { c(v, rep(NA, max.lengthg-length(v)))})
	reptgen=data.frame(do.call(rbind, dupg))
	#write.csv(reptgen,paste("GenotypesRep_",file_name,".csv",sep=""),row.names=F)
}
print("Prepare data...")
#######################################################
## do the filters for missing values in markers
if (missvalG!=0){
	missmark=apply(datos,1,function(y) mean(is.na(y)))
	if(length(which(missmark>missvalG))!=0){
		datos=datos[-which(missmark>missvalG),]
		id=id[-which(missmark>missvalG),]		
	}
}
#######################################################
## do the filters for missing values in genotypes
if (missval!=0){
	missgen=apply(datos,2,function(y) mean(is.na(y)))
	if(length(which(missgen>missval))!=0){datos=datos[,-which(missgen>missval)]}
}else{
	missgen=apply(datos,2,function(y) 2*sum(is.na(y))/(nall*nrow(datos)))
	if(length(which(missgen>0.95))!=0){datos=datos[,-which(missgen>0.95)]}  
}
nacc=ncol(datos)
nalle=nall*nrow(datos)
nmark=nrow(datos)

print("Calculate statistics by individuals...")
#######################################################
datos$pmiss=apply(datos,1,function(y) sum(is.na(y))/nacc)
datos$pmiss1=apply((1-datos),1,function(y) sum(is.na(y))/nacc)
datos$pest=apply(datos[,1:nacc],1,mean,na.rm=T)
datos$pest1=apply((1-datos[,1:nacc]),1,mean,na.rm=T)
id2=data.frame(id,datos$pmiss,datos$pmiss1,datos$pest,datos$pest1)
colnames(id2)=c("Markers", "pmiss", "pmiss1" ,"pest","pest1")
##do the filter for polymorphism
if(mayorque!=0 & menorque!=0){
	exaid=id2[which(id2$pest1<=mayorque & id2$pest<=mayorque & id2$pest>=menorque & id2$pest1>=menorque),]
	datos=datos[which(id2$pest1<=mayorque & id2$pest<=mayorque & id2$pest>=menorque & id2$pest1>=menorque),]
}else{
	exaid=id2
	datos=datos
}
x2=as.data.frame(t(datos[,1:(ncol(datos)-4)]))
idg=as.data.frame(rownames(x2))
colnames(idg)=c("Genotypes")
## with this matrix calculate the genotype he,Ae,ho,etc
tmpmono=which(id2$pest1<=mayorque & id2$pest<=mayorque & id2$pest>=menorque & id2$pest1>=menorque)
#if(length(tmpmono)!=nmark) write.csv(id2[-tmpmono,],paste("MarkOutFilterPoly_",file_name,".csv",sep=""),row.names=F)
######################################################
rm(id2)
rm(tmpmono)
## the ID information useful
id1=exaid[,1]
######################################################
#Calculate specificity 
cond1=which(datos$pest==0|is.na(datos$pest)==TRUE)
cond2=which(datos$pest1==0|is.na(datos$pest1)==TRUE)
longt=apply(datos[,1:nacc],1, length)-apply(datos[,1:nacc],1, function(y) length(which(is.na(y)==TRUE)))
adif01<--(-1)*(apply((datos[,1:nacc]/datos$pest),1, function (y) sum(y*log(y,2),na.rm=T))/longt)
adif02<--(-1)*(apply(((1-datos[,1:nacc])/datos$pest1),1, function (y) sum(y*log(y,2),na.rm=T))/longt)
if (length(cond1)!=0) adif01[cond1]=0
if (length(cond2)!=0) adif02[cond2]=0
##Calculate rareness
rar1=apply(datos[,1:nacc]*adif01,2,sum,na.rm=T)
rar2=apply((1-datos[,1:nacc])*adif02,2,sum,na.rm=T)
rareness<-(rar1+rar2)/nacc
rm(cond1,cond2,longt,rar1,rar2)
######################################################
######################################################
##Begins calculating indexes
######################################################
## observed heterozigosity
nhom=apply(datos[,1:(ncol(datos)-4)], 1, function(y) sum(y==1 | y==0, na.rm=T))
print("Calculate statistics by markers...")
## expected heterozigosity, effective alleles, shannon index, Wright statistics
## prepare new columns for calculate the index PER LOCUS
exaid$pest2=exaid$pest^2
exaid$pest12=exaid$pest1^2
exaid$lpest=log2(exaid$pest)
exaid$lpest1=log2(exaid$pest1)
exaid$shan=-exaid$pest*exaid$lpest
exaid$shan1=-exaid$pest1*exaid$lpest1
exaid$ho=1-nhom/(nacc*(1-(exaid$pmiss)))

## calculus of He, Ho, Ae, shannon index, PER LOCUS
he=1-(exaid$pest2+exaid$pest12)
ae=1/(1-he)
ho=exaid$ho
shannon=exaid$shan+exaid$shan1
inb=he-ho

## To observed distributions indexes (each line correspond and create a plot to each index)
noNA=round(exaid$pmiss,4)
refmark=exaid[,1]
exadiv=data.frame(refmark,he,ho,ae,shannon,adif01,adif02,noNA)
colnames(exadiv)=c("Marker","He","Ho","Ae","Shannon","SpeAllele1","SpeAllele2","%NA")
#write.csv(exadiv,"CalculusPerLocus.csv",row.names=FALSE,quote=FALSE)
rm(exaid)
#############################################################  
He=mean(he)                                										## Expected Heterozygosity(diversidad genetica intrapoblacional)
se_he=sqrt(var(he)/length(he))             										## Standar desviation for Expected Heterozygosity
Ho=mean(ho)                                										## Observed Heterozygosity  (heterocigocidad promedio)
se_ho=sqrt(var(ho)/length(ho))             										## Standar desviation Observed Heterozygosity
Ae=mean(ae)                                										## Number of effective allele (diversidad genetica intrapoblacional)
se_ae=sqrt(var(ae)/length(ae))             										## Standar desviation for Number of effective allele
nm=length(unique(exadiv[,1]))
Shan=sum(shannon)/nm                       										## Shannon diversity index
se_shan=sqrt(var(shannon)/length(shannon)) 										## Standar desviation for Shannon diversity index
Inb=mean(inb)
se_inb=sqrt(var(inb)/length(inb))
PPoli=nm/nmark                                                ## Percentage of polymorphic loci
div=t(t(c(PPoli,He, se_he, Ho, se_ho, Ae, se_ae, Shan, se_shan)))
div=cbind(c("% of polymorphic loci","Expected Heterozygosity","Standar desviation for HE","Observed Heterozygosity", "Standar desviation for HO", 
             "Number of effective allele", "Standar desviation for Ae", "Shannon diversity Index",
             "Standar desviation for ShanIn"),div)

rm(Ae,ae,He,he,ho,Ho,Inb,inb,nhom,PPoli,refmark,se_ae,se_he,se_ho,se_inb,se_shan,Shan,shannon)

nacc1=ncol(x2)
nmark1=nrow(x2)
x2$pmiss=apply(x2,1,function(y) sum(is.na(y))/nacc1)
x2$pmiss1=apply((1-x2),1,function(y) sum(is.na(y))/nacc1)
x2$pest=apply(x2[,1:nacc1],1,mean,na.rm=T)
x2$pest1=apply((1-x2[,1:nacc1]),1,mean,na.rm=T)
id2g=data.frame(idg,x2$pmiss,x2$pmiss1,x2$pest,x2$pest1)
colnames(id2g)=c(colnames(idg), "pmiss", "pmiss1" ,"pest","pest1")
exaidg=id2g
rm(id2g)
## the ID information useful
id1g=exaidg[,1]
## observed heterozigosity
nhomg=apply(x2[,1:nacc1], 1, function(y) sum(y==1 | y==0, na.rm=T))	

## expected heterozigosity, effective alleles, shannon index, Wright statistics
## prepare new columns for calculate the index PER LOCUS
exaidg$pest2=exaidg$pest^2
exaidg$pest12=exaidg$pest1^2
exaidg$lpest=log2(exaidg$pest)
exaidg$lpest1=log2(exaidg$pest1)
exaidg$shan=-exaidg$pest*exaidg$lpest
exaidg$shan1=-exaidg$pest1*exaidg$lpest1
exaidg$ho=1-nhomg/(nacc1*(1-(exaidg$pmiss)))

heg=1-(exaidg$pest2+exaidg$pest12)
aeg=1/(1-heg)
hog=exaidg$ho
shannong=exaidg$shan+exaidg$shan1
## To observed distributions indexes (each line correspond and create a plot to each index)
noNA=round(x2$pmiss,4)
refGenotype=exaidg[,1]
exadivg=data.frame(refGenotype,heg,hog,aeg,shannong,rareness,noNA)
colnames(exadivg)=c("Genotype","He","Ho","Ae","Shannon","rareness","%NA")
#write.csv(exadivg,"CalculusPerGenotype.csv",row.names=FALSE,quote=FALSE)
rm(nacc1,nmark1,x2,exaidg,id1g,nhomg,refGenotype,heg,aeg,hog,shannong,noNA,rareness)

#library(car)
#out="SummaryDiversityAnalysis.csv"
#if("div"%in%ls()==TRUE){
#cat("Diversity","\n","\n",file=out)
#write.table(div, file = out, append = T,quote=F, sep=",",col.names=F,row.names=F)
#}

print("Calculate distance matrix...")
######################################################
##Calculate distances and cluster analysis
######################################################
fr=as.matrix(datos[,1:(ncol(datos)-4)])												## recover the marker information
frn=fr
frn[!is.na(frn)]=1																## no missing values convert to 1
frn[is.na(frn)]=0																  ## missing values convert to 0
N=2*crossprod(frn)																## create square matrix markers information
rm(frn)

if (distk=="Rogers"){  
  #mrdMAT <- mrd_cpp(fr, N)
    aux=matrix(0,nacc,nacc)
	for (i in 1:(nacc-1)) {
		dif <- fr[, i] - fr[, (i+1):nacc, drop = FALSE]
		ss <- colSums(dif^2, na.rm = TRUE)
		aux[i, (i+1):nacc] <- sqrt((2 * ss) / N[i, (i+1):nacc])
	}
	mrdMAT <- aux + t(aux)
}

if (distk=="Nei"){
  fr[is.na(fr)]=0
  aux <- t(fr)%*%fr
  vec <- sqrt(diag(aux))
  aux <- aux/vec[col(aux)]
  aux <- aux/vec[row(aux)]
  aux <- -log(aux)
  mrdMAT=aux
}

colnames(mrdMAT)=colnames(N)
rownames(mrdMAT)=rownames(N)
rm(fr)
gc()

#mrdMAT1=mrdMAT
#colnames(mrdMAT1)=seq(1:dim(mrdMAT)[1])
#DMfs=cbind(ID=seq(1:dim(mrdMAT)[1]),NAME=rownames(mrdMAT),round(mrdMAT1,5))
#write.csv(cbind(ID=seq(1:dim(mrdMAT)[1]),NAME=rownames(mrdMAT),round(mrdMAT1,5)),paste(distk,"Distances.csv",sep=""),quote=FALSE,row.names=FALSE)
#write.table(cbind(NAME=rownames(mrdMAT),round(mrdMAT,5)),"DistancesforGAP.txt",sep="\t",quote=FALSE,row.names=FALSE)
#rm(mrdMAT1)

print("Calculate cluster...")
##graphical representation, the MDS (multidimensional scaling) analysis
mds=cmdscale(mrdMAT, k=3,eig=T)
coord=mds$points
perctCP12=c(round(mds$eig[1]/sum(mds$eig)*100,2),round(mds$eig[2]/sum(mds$eig)*100,2),round(mds$eig[3]/sum(mds$eig)*100,2))
colnames(coord)=c("dim1","dim2","dim3")
#rm(aux,N,mds)
rm(N,mds)

clust=agnes(mrdMAT, method = "ward")

   dend=as.dendrogram(clust)
   coord1=as.data.frame(coord)
   names(coord1)=c("Factor1","Factor2","Factor3")
   #write.csv(coord1,"MDStable.csv",quote=FALSE)
	coord2=cbind(gen=rownames(coord),coord)
	names(coord2)=c("Gen","Factor1","Factor2","Factor3")   

finMark=nrow(datos)
finGen=ncol(datos)-4
   
   ##Create file with the options that were chosen for the analysis
tablaOpt <- data.frame(parametro = character(),valor = character(), stringsAsFactors = FALSE)
# Tipo de datos
if(typedata=="FREQ"){ tablaOpt <- rbind(tablaOpt,data.frame(parametro="InputType", valor="Allele frequencies"))}
if(typedata=="SNP"){ tablaOpt <- rbind(tablaOpt,data.frame(parametro="InputType", valor="SNP"), data.frame(parametro="AA", valor=ht1),
    data.frame(parametro="Aa", valor=ht2), data.frame(parametro="aa", valor=ht3))}
if(typedata=="CUENTA"){ tablaOpt <- rbind(tablaOpt,data.frame(parametro="InputType", valor="Counts"))}
if(typedata=="DistMat"){ tablaOpt <- rbind(tablaOpt,data.frame(parametro="InputType", valor="Distance matrix"))}
# Missing
if(missval!=0){tablaOpt <- rbind(tablaOpt,data.frame(parametro="Missing genotypes (%)", valor=missval*100))}
if(missvalG!=0){tablaOpt <- rbind(tablaOpt,data.frame(parametro="Missing markers (%)", valor=missvalG*100))}
# Otros parámetros
tablaOpt <- rbind(tablaOpt,
  data.frame(parametro="Polymorphism", valor=paste(menorque,"<q<",mayorque)),
  data.frame(parametro="Distance", valor=distk),
  data.frame(parametro="Markers (start)", valor=oriMark),
  data.frame(parametro="Markers (final)", valor=finMark),
  data.frame(parametro="Genotypes (start)", valor=oriGen),
  data.frame(parametro="Genotypes (final)", valor=finGen)
)

	#########################################################################
	#########################################################################
	###For optimize cluster
	#########################################################################
	#########################################################################
	if(gap==TRUE){	
	print("Do optimization...")
		if(methodgap=="gap"){			
			ver=fviz_nbclust(mrdMAT, hcut, k.max=20,method = "gap_stat", nboot = 100)
			gapk1=ver[["data"]]$gap-ver[["data"]]$SE.sim
			test=cbind(ver[["data"]]$gap[-(dim(mrdMAT)[1]-2)],gapk1[-1])
			BestNc=which(test[,1]>=test[,2])  
			if(length(BestNc)>1){ if(min(BestNc)==1){BestNc=BestNc[2]}else{BestNc=min(BestNc)} }
			BestNc=as.numeric(as.character(BestNc))
			if (BestNc==1) {BestNc=1; print("Optimization fail (k=1)")}
			if (is.infinite(BestNc)==T){BestNc=1; print("Optimization fail (k=InF)")}
		}else{
			ver <- fviz_nbclust(mrdMAT,hcut, method = "silhouette",k.max = 20)
			BestNc <- ver$data$clusters[which.max(ver$data$y)]
			if(length(BestNc)>1){ if(min(BestNc)==1){BestNc=BestNc[2]}else{BestNc=min(BestNc)} }
			BestNc=as.numeric(as.character(BestNc))
			if (BestNc==1) {BestNc=1; print("Optimization fail (k=1)")}
			if (is.infinite(BestNc)==T){BestNc=1; print("Optimization fail (k=InF)")}
		}		
	}else{BestNc=nclust}
	
#######################################################
if(mixture==TRUE){
	print("Calculate ancenstry...")
	if (typedata=="SNP" | typedata=="vcfile"){
		genotype_matrix=t(datos[,1:(ncol(datos)-4)])
		rownames(genotype_matrix)=colnames(datos)[1:(ncol(datos)-4)]
		genotype_matrix[genotype_matrix==0]=9
		genotype_matrix[genotype_matrix==0.5]=99
		genotype_matrix[genotype_matrix==1]=999
		genotype_matrix[genotype_matrix==9]=2
		genotype_matrix[genotype_matrix==99]=1
		genotype_matrix[genotype_matrix==999]=0
	}
	if (typedata=="FREQ" | typedata=="CUENTA"){
		genotype_matrix=t(datos[,1:(ncol(datos)-4)])
		rownames(genotype_matrix)=colnames(datos)[1:(ncol(datos)-4)]
		genotype_matrix[genotype_matrix < 0.35]=9
			genotype_matrix[genotype_matrix >= 0.35 & genotype_matrix <= 0.65] <- 99
				genotype_matrix[genotype_matrix > 0.65]=999
				genotype_matrix[genotype_matrix==9]=2
				genotype_matrix[genotype_matrix==99]=1
				genotype_matrix[genotype_matrix==999]=0
		}
		write.geno(genotype_matrix, "datos.geno")
			#########################################################################
			###For do ancestry
			#########################################################################
			#########################################################################	
	if(gap==TRUE & BestNc!=1){
		project <- snmf("datos.geno", K = 1:BestNc, repetitions = 20, entropy=T, project="new",CPU=4)
		# Ver mejor K
		Ks <- 1:BestNc
		ce_list <- lapply(Ks, function(k) {
			cross.entropy(project, K = k)
		})
		ce_matrix <- do.call(cbind, ce_list)
		mean_ce <- colMeans(ce_matrix)
		best_K <- Ks[which.min(mean_ce)]
		best_K=as.numeric(as.character(best_K))
		if(best_K==1){best_K=2}
		cat("Best K is:", best_K, "\n")
		best_run <- which.min(ce_matrix[,best_K])
		cat("Better run for K =", best_K, "is:", best_run, "\n")
		# Extraer matriz Q para K optimo
		qmatrix <- Q(project, K = best_K, run=best_run)	
		rownames(qmatrix) <- rownames(genotype_matrix)
		colnames(qmatrix) <- paste0("Cluster_",1:best_K)
	}else{
		project <- snmf("datos.geno", K = 1:10, repetitions = 20, entropy=T, project="new",CPU=4)
		# Ver mejor K
		Ks <- 1:10
		ce_list <- lapply(Ks, function(k) {
			cross.entropy(project, K = k)
		})
		ce_matrix <- do.call(cbind, ce_list)
		mean_ce <- colMeans(ce_matrix)
		best_K <- Ks[which.min(mean_ce)]
		best_K=as.numeric(as.character(best_K))
		if(best_K==1){best_K=2}		
		cat("Best K is:", best_K, "\n")
		best_run <- which.min(ce_matrix[,best_K])
		cat("Better run for K =", best_K, "is:", best_run, "\n")
		# Extraer matriz Q para K optimo
		qmatrix <- Q(project, K = best_K, run=best_run)	
		rownames(qmatrix) <- rownames(genotype_matrix)
		colnames(qmatrix) <- paste0("Cluster_",1:best_K)
		}
	}else{qmatrix=NULL}

res=list(as.data.frame(div),coord2, getwd(), clust, datos, mrdMAT, perctCP12,BestNc,qmatrix, exadiv, exadivg, tablaOpt)
return(res)
}

