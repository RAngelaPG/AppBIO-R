
rellena <- function(cadena,n=17)
{ 
  tmp <- n-nchar(cadena)
  if(nchar(cadena)<n)  out <- paste(c(cadena,rep(" ",ifelse(nchar(cadena)<7,tmp+6,tmp+2))),collapse="")
  #if(nchar(cadena)<n)  out <- paste(c(cadena,rep(" ",tmp+2)),collapse="")
  if(nchar(cadena)>=n)   out <- paste(substr(cadena,1,n-1),"...",sep="")
  out
}

round2 <- function(numero,digits=0)
{ 
	tmp <- as.character(round(numero,digits)) 
	if(length(grep("\\.",tmp))>0){
		a0 <- unlist(strsplit(tmp,"\\."))[2]
		if(nchar(a0)<digits)	tmp <- paste(c(tmp,rep(0,digits-nchar(a0))),collapse="")
	}else{
		tmp <- paste(c(tmp,".",rep(0,digits)),collapse="")	
	}
	return(as.character(tmp))
}

biplot2 <- function(PC,main="Hola",cex.axis=0.75,cex=0.8)
{	
  a1 <- max(abs(as.vector(PC$loadings[,1])))
  a2 <- max(abs(as.vector(PC$loadings[,2])))
	vari <- PC$sdev^2
	varexpl <- round(100*vari/sum(vari),2)
	xlab <- paste("PC 1 (",varexpl[1],"%)",sep="")
	ylab <- paste("PC 2 (",varexpl[2],"%)",sep="")
	biplot(PC,xlim=c(-a1,a1),cex.axis=cex.axis,cex=cex,xlab=xlab,ylab=ylab,ylim=c(-a2,a2),main=main)# ,xaxt='n',yaxt='n')
	abline(v=0,lty=2)
	abline(h=0,lty=2)	
}

cuadrante2 <- function(punto)
{
		x <- punto[1];y<- punto[2]
		angulo <- abs(atan(y/x)*180/pi)
		limite <- 10
		if(x>=0 & y>=0) out <- ifelse(angulo<limite,4,3)
		if(x<0 & y>=0) out <- ifelse(angulo<limite,2,3)
		if(x<0 & y<0) out <- ifelse(angulo<limite,2,1)
		if(x>=0 & y<0) out <- ifelse(angulo<limite,4,1)
		return(out)
}

biplot3 <- function(PC,main="Hola")
{
	#X <- PC$loadings[,1:2]      ###
	X <- PC[[2]]
	#rownames(X) <- rownames(PC$loadings) ###
	maximo <- apply(X,2,function(x){max(abs(x))})
	for(i in 1:2)  X[,i] <- X[,i]/maximo[i]
	#X <- X/max(abs(as.vector(X)))
	#vari <- PC$sdev^2                     ###
	vari <- PC[[1]]^2
	varexpl <- round(100*vari/sum(vari),2)
	xlab <- paste("PC 1 (",varexpl[1],"%)",sep="")
	ylab <- paste("PC 2 (",varexpl[2],"%)",sep="")
	par(cex=0.9)
	plot(0,0,xlim=c(-1,1)*1.25,ylim=c(-1,1)*1.25,pch="",xlab=xlab,ylab=ylab,main=main)
	par(cex=0.75)
  	for(i in 1:nrow(X)){
		arrows(0,0,X[i,1],X[i,2],code=2,length=0.1,col=2)
		text(X[i,1],X[i,2],rownames(X)[i],col='blue',pos=cuadrante2(X[i,]))
	}
	abline(h=0,lty=2)
	abline(v=0,lty=2)
}

plot_dendrogram <- function(hclust0,main=main,xlab="",ylab="")
{
  par(cex=0.77)
  plot(hclust0,main="",xlab="",hang=-0.15,sub="",yaxt="n",ylab="")
  par(cex=0.9)
  title(xlab=xlab, ylab=ylab, main=main)
}

quita_espacio <- function(cadena)
{
	unlist(lapply(strsplit(cadena," "),function(x)paste(x[1:length(x)],collapse="")))
}

cambia_caracter <- function(cadena)
{
 	especiales <- "-[]=%$+,;<>:/|*?~"
	tmp <- unlist(lapply(strsplit(cadena,""),function(x){
		for(i in 1:length(x))	x[i] <- ifelse(length(grep(x[i],especiales))>0,"_",x[i])
		return(paste(x,collapse=""))
	}))
	return(tmp)
}

putg<-function(cadena){
     noG=which(unlist(lapply(strsplit(cadena,""),function(x){x[1]}))!="g")
     if(length(noG)!=0) cadena=replace(cadena,noG,paste("g",cadena[noG],sep=''))
	 return(cadena)
}

checa_datos <- function(datos,y)
{
	#datos <- datos2 ; y <- traits[j]
	flag <- TRUE
	error <- 0
	mensaje <- "OK"
	datos <- data.frame(datos[,c("Entry","Loc","Rep",y)])
	colnames(datos) <- c("Entry","Loc","Rep",y)
	datos0 <- split(datos,as.character(datos[,"Loc"]))
	calculo <- lapply(datos0,function(x)
	{
		flag <- TRUE
		error <- 0
		mensaje <- "OK"
					
			dd <- split(x,as.character(x[,"Rep"]))
			nRep <- length(dd)
			ee <- lapply(dd,function(xx){
				y0 <- as.numeric(as.character(xx[,y]))
				ID <- as.character(xx[,"Entry"])
				return(data.frame(ID,y0))
			})
			ee0 <- data.frame(ee[[1]])
			if(length(ee)>1) for(ii in 2:length(ee)) ee0 <- merge(ee0,data.frame(ee[[ii]]),by="ID")
			
			# Suma valores no nulos de toda la replica
			sumas <- unlist(lapply(dd,function(xx){
				return(sum(!is.na(as.numeric(as.character(xx[,y])))))
			}))
			SD <- apply(matrix(as.matrix(ee0[,-1]),ncol=ncol(ee0)-1,nrow=nrow(ee0)),2,sd,na.rm=TRUE)
			SD <- SD[!is.na(SD)]
			corre <- cor(matrix(as.matrix(ee0[,-1]),ncol=ncol(ee0)-1,nrow=nrow(ee0)),use="pairwise.complete.obs")
			corre <- corre[lower.tri(corre)]
			corre <- round(corre[!is.na(corre)],3)

			if(length(corre)>0 & any(corre==1) & sum(sumas==0)==0){
				flag <- FALSE
				error <- 1   # correlated
				mensaje <- "Some replicates are correlated"
			}
			if(length(corre)>0 & any(corre==1) & any(sumas==0)){
				flag <- FALSE
				error <- 2   # correlated and null
				mensaje <- "Some replicates are correlated and some others unscored"
			}
			if(any(sumas==0) & length(corre)==0){
				flag <- FALSE
				error <- 3  # null
				mensaje <- "Some replicates are unscored"
			}
			if(sum(sumas)==0 & length(corre)==0){
				flag <- FALSE
				error <- 4  # null
				mensaje <- "All replicates are unscored"
			}
			if(nRep==1){
				flag <- FALSE
				error <- 5  # no replicated
				mensaje <- "There are only one replicate"
			}
			if(length(SD)>0 & sum(SD==0)>1){
				flag <- FALSE
				error <- 6  # valores unicos
				mensaje <- "There are no variance in any replicate"
			}
		return(list(flag=flag,error=error,mensaje=mensaje))
	})

	if(length(calculo)>1)
	{
		tmp <- unlist(lapply(calculo,function(x)x[1]))
		tmp2 <- unlist(lapply(calculo,function(x)x[2]))!=4
		if(sum(tmp2)<2)	## Any with error  4(All records null)
		{
			flag <- FALSE
			mensaje <- paste(sum(!tmp)," out of ",length(tmp)," are unscored locations",sep="")
			error <- 7
		}
		return(list(flag=flag,error=error,mensaje=mensaje))
	}
	if(length(calculo)==1) return(calculo[[1]])
}

################################################
escribe_LOG <- function(linea,FileName="LOG1.txt",directory=dir_root,consola=TRUE)
{
	write.table(paste("   ",paste(linea,collapse=""),sep=""),paste(directory,"/Programs/Rcode/",FileName,sep=""),append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
	if(consola) cat(paste(linea,collapse=""),"\n")
}

elimina_LOG <- function(directory=dir_root,FileName="LOG1.txt")
{
	if(file.exists(paste(directory,"/Programs/Rcode/",FileName,sep="")))
		file.remove(paste(directory,"/Programs/Rcode/",FileName,sep=""))
}

############################
BOXPLOTS <- function()
{
	folderPlot <- "Boxplots"
	if(!file.exists(folderPlot))  dir.create(folderPlot)
	cex0 <- 0.77
	
		index <- which(as.character(datos$Loc)%in%selected)
    		if(length(index)>0) datos <- datos[index,]
  	
  
	for(i in 1:length(traits))
	{
		if("Manag"%in%var.class)
		{
			datos1 <- split(datos,as.character(datos[,"Manag"]))
			for(j in 1:length(datos1))
			{
				tmp <- datos1[[j]]
				tmp <- suppressWarnings(data.frame(Loc=tmp[,"Loc"],y=as.numeric(as.character(tmp[,traits[i]]))))
				to.plot <- split(tmp[,'y'],as.character(tmp[,"Loc"]))
				maxLength <- max(nchar(names(to.plot)))
				flag <- unlist(lapply(to.plot,function(x)any(!is.na(x))))
				if(sum(flag)>0)
				{
				medias <- round(unlist(lapply(to.plot,mean,na.rm=TRUE)),4)
				sdesv <- round(unlist(lapply(to.plot,sd,na.rm=TRUE)),4)
				mins <- suppressWarnings(round(unlist(lapply(to.plot,min,na.rm=TRUE)),4))
				maxs <- suppressWarnings(round(unlist(lapply(to.plot,max,na.rm=TRUE)),4))
				medians <- round(unlist(lapply(to.plot,median,na.rm=TRUE)),4)
				maxs[is.infinite(maxs)] <- NA
				mins[is.infinite(mins)] <- NA
				maximo <- suppressWarnings(max(maxs,na.rm=TRUE))
				minimo <- suppressWarnings(min(mins,na.rm=TRUE))
				ancho <- (maximo-minimo)/5
				mins[is.na(mins)] <- "-"
				maxs[is.na(maxs)] <- "-"
				medias[is.na(medias)] <- "-"
				medians[is.na(medians)] <- "-"
				sdesv[is.na(sdesv)] <- "-"
				
				width <- ifelse(length(to.plot)<5,480,95*length(to.plot))

				if(box_plotFormat=="pdf") pdf(paste(folderPlot,"/Boxplot_",names(datos1)[j],"_",traits[i],".pdf",sep=""),height=7.6,width=7*width/480)
				if(box_plotFormat=="wmf") win.metafile(paste(folderPlot,"/Boxplot_",names(datos1)[j],"_",traits[i],".wmf",sep=""),height=7.6,width=7*width/480)
				if(box_plotFormat=="png") png(paste(folderPlot,"/Boxplot_",names(datos1)[j],"_",traits[i],".png",sep=""),height=520,width=width)

				a <- 2+(1+exp(-0.2*maxLength))*0.45*maxLength
				par(mar=c(a, 4.1, 4.1, 4.1), xpd=FALSE)
				boxplot(to.plot,main=paste("Boxplot. ",names(datos1)[j],sep=""),xlab="",ylab=traits[i],col="orange",
					ylim=c(minimo,maximo+1.30*ancho),yaxt="n",xaxt="n",xlim=c(-0.1,length(to.plot)+0.35))
				abline(h=maximo+0.3*ancho)
				ajuste <- 0.5
				axis(1,1:length(to.plot),labels=names(to.plot),las=2,cex.axis=0.8)
				axis(2,round(seq(minimo,maximo,length=6),1))
				text(c(0.1,1:length(to.plot)),maximo+0.5*ancho,c("StdDesv",sdesv),cex=cex0,adj=c(ajuste,NA))
				text(c(0.1,1:length(to.plot)),maximo+0.73*ancho,c("Max",maxs),cex=cex0,adj=c(ajuste,NA))
				text(c(0.1,1:length(to.plot)),maximo+0.96*ancho,c("Median",medians),cex=cex0,adj=c(ajuste,NA))
				text(c(0.1,1:length(to.plot)),maximo+1.19*ancho,c("Mean",medias),cex=cex0,adj=c(ajuste,NA))
				text(c(0.1,1:length(to.plot)),maximo+1.42*ancho,c("Min",mins),cex=cex0,adj=c(ajuste,NA))
				for(k in 1:length(to.plot))	rect(k-0.5,maximo+0.3*ancho,k-0.5,maximo+1.55*ancho)
				dev.off()
				}
			}
		}

		tmp <- suppressWarnings(data.frame(Loc=datos[,"Loc"],y=as.numeric(as.character(datos[,traits[i]]))))
		to.plot <- split(tmp[,'y'],as.character(tmp[,"Loc"]))
		maxLength <- max(nchar(names(to.plot)))
		flag <- unlist(lapply(to.plot,function(x)any(!is.na(x))))
		if(sum(flag)>0)
		{
			medias <- round(unlist(lapply(to.plot,mean,na.rm=TRUE)),4)
			sdesv <- round(unlist(lapply(to.plot,sd,na.rm=TRUE)),4)
			mins <- suppressWarnings(round(unlist(lapply(to.plot,min,na.rm=TRUE)),4))
			maxs <- suppressWarnings(round(unlist(lapply(to.plot,max,na.rm=TRUE)),4))
			medians <- round(unlist(lapply(to.plot,median,na.rm=TRUE)),4)
			maximo <- max(maxs,na.rm=TRUE)
			minimo <- min(mins,na.rm=TRUE)
			ancho <- (maximo-minimo)/5
			maxs[is.infinite(maxs)] <- NA
			mins[is.infinite(mins)] <- NA
			mins[is.na(mins)] <- "-"
			maxs[is.na(maxs)] <- "-"
			medias[is.na(medias)] <- "-"
			medians[is.na(medians)] <- "-"
			sdesv[is.na(sdesv)] <- "-"

			width <- ifelse(length(to.plot)<5,480,95*length(to.plot))
			if(box_plotFormat=="pdf") pdf(paste(folderPlot,"/Boxplot_AllLoc_",traits[i],".pdf",sep=""),height=7.6,width=7*width/480)
			if(box_plotFormat=="wmf") win.metafile(paste(folderPlot,"/Boxplot_AllLoc_",traits[i],".wmf",sep=""),height=7.6,width=7*width/480)
			if(box_plotFormat=="png") png(paste(folderPlot,"/Boxplot_AllLoc_",traits[i],".png",sep=""),height=520,width=width)

			a <- 2+(1+exp(-0.2*maxLength))*0.45*maxLength
			par(mar=c(a, 4.1, 4.1, 4.1), xpd=FALSE)
			boxplot(to.plot,main="Boxplot. All Locations",xlab="",ylab=traits[i],col="orange",
			ylim=c(minimo,maximo+1.30*ancho),yaxt="n",xaxt="n",xlim=c(-0.1,length(to.plot)+0.35))
			abline(h=maximo+0.3*ancho)
			ajuste <- 0.5
			axis(1,1:length(to.plot),labels=names(to.plot),las=2,cex.axis=0.8)
			axis(2,round(seq(minimo,maximo,length=6),1))
			text(c(0.1,1:length(to.plot)),maximo+0.5*ancho,c("StdDesv",sdesv),cex=cex0,adj=c(ajuste,NA))
			text(c(0.1,1:length(to.plot)),maximo+0.73*ancho,c("Max",maxs),cex=cex0,adj=c(ajuste,NA))
			text(c(0.1,1:length(to.plot)),maximo+0.96*ancho,c("Median",medians),cex=cex0,adj=c(ajuste,NA))
			text(c(0.1,1:length(to.plot)),maximo+1.19*ancho,c("Mean",medias),cex=cex0,adj=c(ajuste,NA))
			text(c(0.1,1:length(to.plot)),maximo+1.42*ancho,c("Min",mins),cex=cex0,adj=c(ajuste,NA))
			for(k in 1:length(to.plot))	rect(k-0.5,maximo+0.3*ancho,k-0.5,maximo+1.55*ancho)
			dev.off()
		}
		
	}
	escribe_LOG(" ")
	escribe_LOG("Boxplots done")
}

############################################################
BASIC_STATS <- function()
{
	out <- c()
	for(i in 1:length(traits))
	{
		if("Manag"%in%var.class)
		{
			datos1 <- split(datos,as.character(datos[,"Manag"]))
			for(j in 1:length(datos1))
			{
				tmp <- datos1[[j]]
				tmp <- suppressWarnings(data.frame(Loc=tmp[,"Loc"],y=as.numeric(as.character(tmp[,traits[i]]))))
				to.plot <- split(tmp[,'y'],as.character(tmp[,"Loc"]))
				medias <- round(unlist(lapply(to.plot,mean,na.rm=TRUE)),4)
				conteo <- round(unlist(lapply(to.plot,length)),4)
				nNA <- round(unlist(lapply(to.plot,function(x)sum(is.na(x)))),4)
				sdesv <- round(unlist(lapply(to.plot,sd,na.rm=TRUE)),4)
				mins <- suppressWarnings(round(unlist(lapply(to.plot,min,na.rm=TRUE)),4))
				maxs <- suppressWarnings(round(unlist(lapply(to.plot,max,na.rm=TRUE)),4))
				quantile25 <- round(unlist(lapply(to.plot,quantile,prob=0.25,na.rm=TRUE)),4)
				medians <- round(unlist(lapply(to.plot,median,na.rm=TRUE)),4)
				quantile75 <- round(unlist(lapply(to.plot,quantile,prob=0.75,na.rm=TRUE)),4)
				tmp <- cbind(conteo,nNA,mins,quantile25,medians,quantile75,maxs,medias,sdesv)
				tmp[is.infinite(tmp) | is.na(tmp)] <- "-"
				tmp <- cbind(rep(traits[i],length(to.plot)),rep(names(datos1)[j],length(to.plot)),names(mins),tmp)
				out <- rbind(out,tmp)
			}
		}
		tmp <- suppressWarnings(data.frame(Loc=datos[,"Loc"],y=as.numeric(as.character(datos[,traits[i]]))))
		to.plot <- split(tmp[,'y'],as.character(tmp[,"Loc"]))
		medias <- round(unlist(lapply(to.plot,mean,na.rm=TRUE)),4)
		conteo <- round(unlist(lapply(to.plot,length)),4)
		nNA <- round(unlist(lapply(to.plot,function(x)sum(is.na(x)))),4)
		sdesv <- round(unlist(lapply(to.plot,sd,na.rm=TRUE)),4)
		mins <- suppressWarnings(round(unlist(lapply(to.plot,min,na.rm=TRUE)),4))
		maxs <- suppressWarnings(round(unlist(lapply(to.plot,max,na.rm=TRUE)),4))
		quantile25 <- round(unlist(lapply(to.plot,quantile,prob=0.25,na.rm=TRUE)),4)
		medians <- round(unlist(lapply(to.plot,median,na.rm=TRUE)),4)
		quantile75 <- round(unlist(lapply(to.plot,quantile,prob=0.75,na.rm=TRUE)),4)
		tmp <- cbind(conteo,nNA,mins,quantile25,medians,quantile75,maxs,medias,sdesv)
		tmp[is.infinite(tmp) | is.na(tmp)] <- "-"
		tmp <- cbind(rep(traits[i],length(to.plot)),rep("-",length(to.plot)),names(mins),tmp)
		out <- rbind(out,tmp,rep("",ncol(tmp)))
	}
	colnames(out) <- c("Trait","Management","Location","N","N_Missing","Minimum","Quantile25","Median","Quantile75",
	"Maximum","Mean","Std Desv")
	write.table(out,"Basic_Statistics.csv",sep=",",quote=FALSE,row.names=FALSE)
	escribe_LOG(" ")
	escribe_LOG("Summary calculated")
}

#########################################################
HISTOGRAMAS <- function()
{
	folderPlot <- "Histograms"
	if(!file.exists(folderPlot))  dir.create(folderPlot)
	cex0 <- 0.77
	for(i in 1:length(traits))
	{
		if("Manag"%in%var.class)
		{
			datos1 <- split(datos,as.character(datos[,"Manag"]))
			for(j in 1:length(datos1))
			{
				tmp <- datos1[[j]]
				tmp <- suppressWarnings(data.frame(Loc=tmp[,"Loc"],y=as.numeric(as.character(tmp[,traits[i]]))))
				to.plot <- split(tmp[,'y'],as.character(tmp[,"Loc"]))
				flag <- unlist(lapply(to.plot,function(x)any(!is.na(x))))
				if(sum(flag)>0)
				{
					for(k in 1:length(to.plot))
					{
						main <- paste("Histogram. ",names(datos1)[j],"\n",names(to.plot)[k],sep="")
						if(flag[k])
						{	
							if(histogramasFormat=="pdf") pdf(paste(folderPlot,"/Histogram_",names(datos1)[j],"_",names(to.plot)[k],"_",traits[i],".pdf",sep=""))
							if(histogramasFormat=="wmf") win.metafile(paste(folderPlot,"/Histogram_",names(datos1)[j],"_",names(to.plot)[k],"_",traits[i],".wmf",sep=""))
							if(histogramasFormat=="png") png(paste(folderPlot,"/Histogram_",names(datos1)[j],"_",names(to.plot)[k],"_",traits[i],".png",sep=""))

							hist(to.plot[[k]],main=main,xlab=traits[i],col="orange")
							dev.off()
						}
					}
				}
			}
		}

		tmp <- suppressWarnings(data.frame(Loc=datos[,"Loc"],y=as.numeric(as.character(datos[,traits[i]]))))
		to.plot <- split(tmp[,'y'],as.character(tmp[,"Loc"]))
		flag <- unlist(lapply(to.plot,function(x)any(!is.na(x))))
		if(sum(flag)>0)
		{
			for(k in 1:length(to.plot))
			{
				if(flag[k]){
					if(histogramasFormat=="pdf") pdf(paste(folderPlot,"/Histogram_",names(to.plot)[k],"_",traits[i],".pdf",sep=""))
					if(histogramasFormat=="wmf") win.metafile(paste(folderPlot,"/Histogram_",names(to.plot)[k],"_",traits[i],".wmf",sep=""))
					if(histogramasFormat=="png") png(paste(folderPlot,"/Histogram_",names(to.plot)[k],"_",traits[i],".png",sep=""))

					hist(to.plot[[k]],main=paste("Histogram. ",names(to.plot)[k],sep=""),xlab=traits[i],col="orange")
					dev.off()
				}
			}
		}
	}
	
	escribe_LOG(" ")
	escribe_LOG("Histograms done")
}

#==============================

get_h2 <- function(datos=datos2,traits=traits,whichTrait=i,covariate=covariate)
{      
	envirs <- unique(as.character(datos[,"Loc"]))
	h2 <- rep(NA,length(envirs))
	names(h2) <- envirs

	datos.tmp <- split(datos,as.character(datos$Loc))
	tmp <- lapply(datos.tmp,function(x){
		for(k in 1:length(var.class))  x[,var.class[k]] <- as.factor(as.character(x[,var.class[k]]))	
		out <- rep(NA,5)
		checa <- suppressWarnings(checa_datos(x,traits[whichTrait]))
		if(checa$flag)
		{
			tt <- suppressWarnings(data.frame(x[,var.class],y=as.numeric(as.character(x[,traits[whichTrait]])))) 
			tt <- data.frame(tt,y_est=scale(tt$y))
			flagCOV <- FALSE
			if(!is.null(covariate))
			{	
				COV <- suppressWarnings(apply(matrix(as.matrix(x)[,covariate],ncol=length(covariate)),2,as.numeric))
				COV <- scale(COV)
				colnames(COV) <- paste("Cov",1:length(covariate),sep="")
				checaCOV <- list()
				for(k in 1:length(covariate)) checaCOV[[k]] <- suppressWarnings(checa_datos(x,covariate[k]))
				flagCOV2 <- unlist(lapply(checaCOV,function(x)x$flag))
				flagCOV <- prod(flagCOV2)==1
				if(flagCOV) tt <- data.frame(tt,COV)
			}

			fm.text <- ifelse(ExpDes=="RCB","lmer(y~(1|Entry)+(1|Rep)","lmer(y~(1|Entry)+(1|Rep)+(1|Block:Rep)")						
			tmpCOV <- ifelse(flagCOV & !is.null(covariate) & whichTrait==1,paste("+Cov",1:length(covariate),collapse="",sep=""),"")
			fm.text <- paste(fm.text,tmpCOV,",data=tt)",sep="")
			fm <- suppressWarnings(eval(parse(text=fm.text)))

			#  Compute h2
			varcorr <- VarCorr(fm)
			varG <- as.vector(varcorr$'Entry')
			varErr <- attr(varcorr,'sc')^2
			nRep <- length(unique(as.character(x$Rep)))
			h2 <- round(varG/(varG + varErr/nRep),3)
			out[1:5] <- c(as.character(unique(x$Loc)),round(varG,3),round(varErr,3),h2,"")
		}
		if(!checa$flag) out[c(1,5)] <- c(as.character(unique(x$Loc)),checa$mensaje)
		return(out)		
	})		
	out <- do.call('rbind',tmp)
	colnames(out) <- c("Location","Genotype Variance","Residual Variance","Heritability","")	
	return(out)	
}



