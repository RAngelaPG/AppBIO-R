server <- function(input,output,session) {
  
  if (!interactive()) {
    session$onSessionEnded(function() {
      stopApp()
      q("no")
    })
  }
  
  dirappViT=getwd()
  dirappViT=str_replace(dirappViT,"shiny","")
  print(dirappViT)	   

output$menu <- renderMenu({
    sidebarMenu(
      menuItem(tagList(span("Data",title="Load files")),tabName="datos", icon = icon("database") 
      ),      
      menuItem(tagList(img(src="Biopng.png",height="30"),
                      span("Biodiversity analysis", title="Do the biodiversity analysis, in order to calculate 
heterozygosity, diversity among and within groups, shannon index, number 
of effective allele, percent of polymorphic loci, Rogers distance, Nei 
distance, cluster analysis and multidimensional scaling 2D plot and 3D 
plot; you can included external groups for colored the dendogram or MDS 
plots.")),tabName="bio"),			   
	  menuItem(tagList(img(src="Biopng.png",height="30"),
                      span("CoreSubset", title="Obtain sampling core subsets from genetic resources while maintaining
as much as possible the genetic diversity of the original collection")),tabName="core"),
	  menuItem(tagList(img(src="Biopng.png",height="30"),
                      span("GWAS", title="Research method that compares the DNA of a large group of individuals 
					  to identify genetic variations associated with a disease or trait")),tabName="gwas")			   			   
               #menuSubItem("Biodivesity analysis", tabName="bio", icon =  icon("envira")),
               #menuSubItem("CoreSubset", tabName="core", icon = icon("envira")),
			   #menuSubItem("GWAS", tabName="gwas", icon = icon("envira"))
      )
    #)
    
  })
  
##################################################################################################################################
  #Genetic data#
  ##################################################################################################################################
  shinyFileChoose(input, 'filegen', roots = getVolumes(),filetypes=c('', 'csv','vcf','gz'))
  shinyFileChoose(input, 'fileRdata', roots = getVolumes(),filetypes=c('', 'RData'))
  
  observe({	
    if(input$typedata=="SNP"){
      shinyjs::show(id="ht1")
      shinyjs::show(id="ht2")
      shinyjs::show(id="ht3")
    }else{
      shinyjs::hide(id="ht1")
      shinyjs::hide(id="ht2")
      shinyjs::hide(id="ht3")
    }
	
	if(input$startAna=="StarChrom"){
		shinyjs::hide(id="typedata")
		shinyjs::hide(id="missval")
		shinyjs::hide(id="mayorque")
		shinyjs::hide(id="menorque")	
		shinyjs::show(id="fileVCF")
		shinyjs::hide(id="ht1")
		shinyjs::hide(id="ht2")
		shinyjs::hide(id="ht3")
		updateRadioButtons(session, "typedata", selected = "vcfile")		
	}	
	if(input$startAna=="pavs"){
		shinyjs::hide(id="typedata")
		shinyjs::hide(id="missval")
		shinyjs::hide(id="mayorque")
		shinyjs::hide(id="menorque")	
		shinyjs::hide(id="fileVCF")
		shinyjs::hide(id="ht1")
		shinyjs::hide(id="ht2")
		shinyjs::hide(id="ht3")		
	}
	if(input$startAna=="StarBio"){
		shinyjs::show(id="typedata")
		shinyjs::show(id="missval")
		shinyjs::show(id="mayorque")
		shinyjs::show(id="menorque")
		shinyjs::hide(id="fileVCF")		
	}	
  })
  
  
	GenInfo <- reactiveValues(dfgen=NULL,posit=NULL,UploadRd=NULL,hapmap=NULL)
	observeEvent(input$readgenofile, {
  #myDataGen <- reactive({
    #SelFile="Data"
	UploadRd=NULL
	#Gdata=0
	inFilegen=parseFilePaths(roots=getVolumes(), input$filegen)
	inFileRdata=parseFilePaths(roots=getVolumes(), input$fileRdata)	
	 if(nrow(inFilegen)!=0 & nrow(inFileRdata)==0){Gdata=1;SelFile="Data"}
	 if(nrow(inFilegen)==0 & nrow(inFileRdata)!=0){Gdata=1;SelFile="RData"}
	 if(nrow(inFilegen)==0 & nrow(inFileRdata)==0){Gdata=0}
	 if(nrow(inFilegen)!=0 & nrow(inFileRdata)!=0){NAdata=0}else{NAdata=1} 
	 
	 #SNP-17
	file_attr_SNP=c("AlleleID","CloneID","AlleleSequence","SNP","SnpPosition","CallRate","OneRatioRef","OneRatioSnp",
                "FreqHomRef","FreqHomSnp","FreqHets","PICRef","PICSnp","AvgPIC","AvgCountRef","AvgCountSnp","RepAvg")
	#COUNTS-32
	file_attr_COUNTS=c("AlleleID","CloneID","ClusterTempIndex","AlleleSequence","ClusterConsensusSequence","ClusterSize",
                   "AlleleSeqDist","SNP","SnpPosition","CallRate","OneRatioRef","OneRatioSnp","FreqHomRef","FreqHomSnp",
                   "FreqHets","PICRef","PICSnp","AvgPIC","AvgCountRef","AvgCountSnp","RatioAvgCountRefAvgCountSnp",
                   "FreqHetsMinusFreqMinHom","AlleleCountsCorrelation","aggregateTagsTotal","DerivedCorrMinusSeedCorr",
                   "RepRef","RepSNP","RepAvg","PicRepRef","PicRepSNP","TotalPicRepRefTest","TotalPicRepSnpTest")
	#PAVS-15
	file_attr_PAVS=c("CloneID","ClusterTempIndex","AlleleSequence","ClusterConsensusSequence","ClusterSize","CallRate","OneRatio",
               "PIC","AvgReadDepth","StDevReadDepth","Qpmr","aggregateTagsTotal","Reproducibility","PicRep","TotalPicRepTest")
    validate(
      need(Gdata != 0, "Please upload data"),
	  need(NAdata != 0, "Both data selected, please close and try again")	  
    )
	withProgress(message = 'Getting...', value = 0,{
	incProgress(1/2, detail = "Wait, Please!")
		
		if(SelFile=="Data"){
			if(input$startAna=="pavs"){
				dfgen = fread(as.character(inFilegen$datapath),header = TRUE,sep = ",",na=c("NA",".","-",""))
				posit=NULL
				hapmap=NULL
				if(length(intersect(file_attr_PAVS, dfgen[4,])) == 15){
					gename=as.matrix(dfgen[4,16:dim(dfgen)[2]])
					idname=as.character(as.matrix(dfgen[-c(1:4),1]))
					dfgen=as.data.frame(dfgen[-c(1:4),-c(1:15)])
					dfgen=apply(dfgen,2,as.numeric)
					colnames(dfgen)=gename
					rownames(dfgen)=idname
				}else{
					gename=names(dfgen)[-1]
					idname=as.character(as.matrix(dfgen[,1]))
					dfgen=as.data.frame(dfgen[,-1])
					dfgen=apply(dfgen,2,as.numeric)
					colnames(dfgen)=gename
					rownames(dfgen)=idname
				}
			}else{
			if(input$typedata=="vcfile"){
				filevcf=read_vcf(as.character(inFilegen$datapath), ploidity = 2, na_reps = c("./.","-","NN","NA"))
				hapmap=filevcf[[2]]
				filevcf=filevcf[[1]]				
				posit=data.frame(cbind(CHROM=as.character(filevcf@chromosome),POS=filevcf@position,ID=filevcf@loc.names))		
				dfgen=t(as.matrix(filevcf))
				class(dfgen) = "numeric"
				dfgen=data.frame(cbind(rs=rownames(dfgen),data.frame(dfgen)))				
			}else{
			if(input$typedata=="CUENTA"){
				dfgen = read.csv.ffdf(file=as.character(inFilegen$datapath),VERBOSE=T)
				posit=NULL
				hapmap=NULL
				if(length(intersect(file_attr_COUNTS, as.matrix(dfgen[4,]))) == 32){
					gename=c("AlleleID",as.matrix(dfgen[4,33:dim(dfgen)[2]]))
					coluno=dfgen[-c(1:4),1]					
					dfgen=as.data.frame(dfgen[-c(1:4),-c(2:32)])
					dfgen <- apply(dfgen,2,as.numeric)
					dfgen <-cbind(coluno,dfgen)					
					colnames(dfgen)=gename
					dfgen=as.ffdf(dfgen)
				}
			}else{
				dfgen = fread(as.character(inFilegen$datapath),header = TRUE,sep = ",",na=c("NA",".","-",""))
				posit=NULL
				hapmap=NULL
				if(input$typedata=="SNP"){
				if(length(intersect(file_attr_SNP, dfgen[4,])) == 17){
					gename=c("AlleleID",as.matrix(dfgen[4,18:dim(dfgen)[2]]))
					coluno=dfgen[-c(1:4),1]
					dfgen <-as.data.frame(dfgen[-c(1:4),-c(1:17)])
					dfgen <- apply(dfgen,2,as.numeric)
					dfgen <-cbind(coluno,dfgen)
					colnames(dfgen)=gename
				}
			  }
			}
			}
			
			}
		}else{
			dfgen=NULL
			posit=NULL
			hapmap=NULL
			load(as.character(inFileRdata$datapath),UploadRd<-new.env())
			vars=names(UploadRd$Aux[[1]])			
			#Actualiza el selectinput de acuerdo a la base de datos cargada
			updateSelectInput(session,'xcol', 'X Variable',choices = vars,selected=vars[2])
			updateSelectInput(session,'ycol', 'Y Variable',choices = vars,selected=vars[3])
			updateSelectInput(session,'zcol', 'Z Variable',choices = vars,selected=vars[4])
			updateSelectInput(session,'catv', 'Group',choices = vars,selected=vars[5])
			updateSelectInput(session,'eti', 'Label',choices = vars,selected=vars[1])
			updateSelectInput(session,'xcol3D', 'X Variable',choices = vars,selected=vars[2])
			updateSelectInput(session,'ycol3D', 'Y Variable',choices = vars,selected=vars[3])
			updateSelectInput(session,'zcol3D', 'Z Variable',choices = vars,selected=vars[4])
			updateSelectInput(session,'catv3D', 'Group',choices = vars,selected=vars[5])
			updateSelectInput(session,'eti3D', 'Label',choices = vars,selected=vars[1])
			updateSelectInput(session,'catvdend', 'Group',choices = vars,selected=vars[5])    
			updateTextInput(session,'tx','X Axis Label',value = paste0('Factor 1 (',UploadRd$DivAna[[7]][1],'%)'))
			updateTextInput(session,'ty','Y Axis Label',value = paste0('Factor 2 (',UploadRd$DivAna[[7]][2],'%)'))
			updateTextInput(session,'tz','Z Axis Label',value = paste0('Factor 3 (',UploadRd$DivAna[[7]][3],'%)'))
			updateTextInput(session,'tx3D','X Axis Label',value = paste0('Factor 1 (',UploadRd$DivAna[[7]][1],'%)'))
			updateTextInput(session,'ty3D','Y Axis Label',value = paste0('Factor 2 (',UploadRd$DivAna[[7]][2],'%)'))
			updateTextInput(session,'tz3D','Z Axis Label',value = paste0('Factor 3 (',UploadRd$DivAna[[7]][3],'%)'))
			updateTextInput(session,'nclust','No. Clusters',value = paste0(UploadRd$DivAna[[8]]))
		}
    incProgress(1, detail = "Finish")
	Sys.sleep(1)
	})
	
	if(SelFile=="RData"){
		shinyjs::hide(id="fileenvbio")
		shinyjs::hide(id="quitomono")
		shinyjs::hide(id="gapS")
		shinyjs::hide(id="nclust")
		shinyjs::hide(id="typedata")
		shinyjs::hide(id="distk")
	}else{
		shinyjs::show(id="fileenvbio")
		shinyjs::show(id="quitomono")
		shinyjs::show(id="gapS")
		shinyjs::show(id="nclust")
		shinyjs::show(id="distk")
	}
	
    GenInfo$dfgen<-dfgen
	GenInfo$posit<-posit
	GenInfo$UploadRd<-UploadRd
	GenInfo$hapmap<-hapmap
	
	#return(list(dfgen,posit,SelFile,UploadRd))
  })
  #Ver datos en tabla dinamica
  output$seeDataGen<-DT::renderDataTable({	
    if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(!is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){NAdata=0}else{NAdata=1} 
	 if(input$startAna=="StarBio"){go1=1}else{go1=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(NAdata != 0, "Both data selected, please close and try again"),	  
	  need(go1 != 0, "go to the correspond tab of type analysis")		  
    )
	if(input$startAna=="StarBio"){		
		if(SelFile=="RData"){
			#load(as.character(parseFilePaths(roots=getVolumes(), input$fileRdata)$datapath))
			UpRD=GenInfo$UploadRd
			seedatos<-UpRD$DivAna[[5]][,1:(ncol(UpRD$DivAna[[5]])-4)]
		}else{
			seedatos<-as.data.frame(GenInfo$dfgen)
		}
	}
    datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 15,width="100%", scrollX = TRUE))
	
  })
  
##################################################################################################################################################################  
  #See heatmap PAVs
 ##################################################################################################################################################################
  output$seePAVS<-renderPlotly({
  
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(!is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){NAdata=0}else{NAdata=1} 
	 if(input$startAna=="pavs"){go1=1}else{go1=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(NAdata != 0, "Both data selected, please close and try again"),	  
	  need(go1 != 0, "go to the correspond tab of type analysis")		  
    )
	if(input$startAna=="pavs"){
	#withProgress(message = 'Getting...', value = 0,{
	#	incProgress(1/2, detail = "Wait, Please!")		
			
		if(SelFile=="Data"){
			inFilegen=parseFilePaths(roots=getVolumes(), input$filegen)
			dirfile=as.character(inFilegen$datapath)
			filename=as.character(inFilegen$name)
			outFolder <- cambia_caracter(paste0("Output_PAVs_",filename))
			filename=str_replace(dirfile,filename,"")			
			setwd(filename)
			if(!file.exists(outFolder)) dir.create(outFolder)
			setwd(outFolder)			
			use=GenInfo$dfgen
			save(use,file="PAVS.RData")	
			saveurl=getwd()
			pavsplot=plot_ly(x=colnames(use)[1:dim(use)[2]],y=rownames(use)[1:dim(use)[1]],z = use[1:dim(use)[1],1:dim(use)[2]], colorscale="Jet",type = "heatmap")%>%
			layout(xaxis = list(showticklabels = T), yaxis = list(showticklabels = F))
			withr::with_dir(file.path(saveurl),saveWidget(pavsplot,'PAVSPlot.html', selfcontained = F))
			pavsplot
		}else{
			#use=myDataGen()[[4]]
			inFileRdata=parseFilePaths(roots=getVolumes(), input$fileRdata)	
			load(as.character(inFileRdata$datapath))
			pavsplot=plot_ly(x=colnames(use)[1:dim(use)[2]],y=rownames(use)[1:dim(use)[1]],z = use[1:dim(use)[1],1:dim(use)[2]], colorscale="Jet",type = "heatmap")%>%
			layout(xaxis = list(showticklabels = T), yaxis = list(showticklabels = F))			
			pavsplot
		}
	#	incProgress(1, detail = "Finish")
	#	Sys.sleep(1)
	#	})
	}else{
			fig = plotly::plot_ly()
            fig = fig %>% plotly::add_annotations(text = "go to the correspond tab of type analysis", x = 1, y = 1)#
            fig
	}
	
  })
 
##################################################################################################################################################################
  #Ver grafico de cromosomas
##################################################################################################################################################################
  output$seeChromPlot<-renderChromoMap({
	inFilegen=parseFilePaths(roots=getVolumes(), input$filegen)
	inFileRdata=parseFilePaths(roots=getVolumes(), input$fileRdata)	
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(!is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){NAdata=0}else{NAdata=1}
	 if(input$typedata=="vcfile"){ dataT1=1 }else {dataT1=0}	 
	 if(input$startAna=="StarChrom"){go1=1}else{go1=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(NAdata != 0, "Both data selected, please close and try again"),
	  need(dataT1 != 0, "Option only available for VCF files"),		  	  	  
	  need(go1 != 0, "go to the correspond tab of type analysis")		  
    )
	
	
	if(input$startAna=="StarChrom"){	
		withProgress(message = 'Getting...', value = 0,{
		incProgress(1/2, detail = "Wait, Please!")		
		if(SelFile=="Data"){
			dirfile=as.character(inFilegen$datapath)    		
				#if(input$startAna=="StarBio"){
					posit=data.frame(GenInfo$posit)
				#}else{
				#	posit=vcfR::read.vcfR(as.character(dirfile))
				#	posit=data.frame(posit@fix)
				#}
				tmark=dim(posit)[1]
				posit$POS=as.numeric(posit$POS)
				posit$CHROM=as.character(posit$CHROM)
				posit$ID=as.character(posit$ID)
				tabamv1=data.frame(table(as.factor(posit$CHROM)))				
				colnames(tabamv1)=c("Chromosome","NMarks")
				dimt1=dim(tabamv1)[1]
				#chromi=as.character(tabamv1[which(tabamv1$NMarks>=(sum(tabamv1$NMarks)*0.05)),1])
				chromi=as.character(tabamv1[which(tabamv1$NMarks>=1),1])				
				chromi1=chromi[grep("chr",chromi)]
				if(length(chromi1)==0){chromi=chromi}else{chromi=chromi1}
				if(length(which(posit$CHROM%in%chromi==F))!=0){posit$CHROM[which(posit$CHROM%in%chromi==F)]=NA}		
				if(length(which(is.na(posit[,1])))!=0){posit=posit[-which(is.na(posit[,1])),]}
				tabamv1=data.frame(table(as.factor(posit$CHROM)))		
				colnames(tabamv1)=c("Chromosome","NMarks")
				if(dimt1>dim(tabamv1)[1]) {
				shinyalert(title = "Important message", 
					text="Extra chromosome information was \n found in the blast and was deleted",closeOnEsc = FALSE, 
					type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
				)
					#print("Extra chromosome information was found in the blast and was deleted") 
				}
				filename=as.character(inFilegen$name)
				filenamef=str_replace(filename,".csv","")
				filename=strsplit(filename,"\\.")[[1]][1]
				out=paste0("SummaryBlast_",filename,".csv")
				
				outFolder <- cambia_caracter(paste0("ChromMap_",filename))
				setwd(str_replace(dirfile,filenamef,""))
				if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
				setwd("Output_BIO-R")
				if(!file.exists(outFolder)) dir.create(outFolder)
				setwd(outFolder)			
			
				cat("\n","Total markers: ",tmark,"\n","\n",file=out,append=T)
				cat("\n","Markers in blast: ",dim(posit)[1],"\n","\n",file=out,append=T)
				cat("\n","Markers in blast for each chromosome: ","\n","\n",file=out,append=T)
				write.table(tabamv1, file = out, append = T,quote=F, sep=",",col.names=T,row.names=F)
				system2('open',args=out,wait=F)
					
				chromosome_file=aggregate(POS ~ CHROM, data = posit, max)
				chromosome_file=data.frame(cbind(CHROM=chromosome_file$CHROM,STAR=1,END=chromosome_file$POS))
				chromosome_file$STAR=as.numeric(chromosome_file$STAR)
				chromosome_file$END=as.numeric(chromosome_file$END)
		
				annotation_file=data.frame(cbind(ID=posit[,"ID"],CHROM=posit[,"CHROM"],STAR=posit[,"POS"],END=posit[,"POS"]+10))
				annotation_file$STAR=as.numeric(annotation_file$STAR)
				annotation_file$END=as.numeric(annotation_file$END)
				
				chrdim=split(annotation_file,annotation_file$CHROM)
				minchr=lapply(chrdim,function(x){min(min(x$STAR),min(x$END))})
				maxchr=lapply(chrdim,function(x){max(max(x$STAR),max(x$END))})
				chrdim=data.frame(cbind(minchr=unlist(minchr),maxchr=unlist(maxchr)))
				write.csv(chrdim,"ChromDim.csv")
		
				save(chromosome_file,annotation_file,file="MapChrom.RData")	
			}else{
				load(as.character(inFileRdata$datapath))
			}			
		incProgress(1, detail = "Finish")
		Sys.sleep(1)
		})	
	}
	
	chromoMap(list(chromosome_file),list(annotation_file),n_win.factor = 5, export.options=T, fixed.window=F, remove.last.window=T,plot.shift=T, id="chrom1")
	
  })
  
  shinyFileChoose(input, 'fileVCF', roots = getVolumes(),filetypes=c('vcf','gz'))

##################################################################################################################################################################
  #Ver grafico comparacion de blast cromosomas
##################################################################################################################################################################
  output$BChromPlot<-renderChromoMap({
	 inFileRdata=parseFilePaths(roots=getVolumes(), input$fileRdata)	
	 if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(!is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){NAdata=0}else{NAdata=1}
	 if(input$startAna=="StarChrom"){go1=1}else{go1=0}
	 
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(NAdata != 0, "Both data selected, please close and try again"),	  
	  need(go1 != 0, "go to the correspond tab of type analysis")		  
    )
	 
	 inFileRdata=parseFilePaths(roots=getVolumes(), input$fileRdata)	
	 inFilegen2=parseFilePaths(roots=getVolumes(), input$fileVCF)
	 if(nrow(inFilegen2)!=0 & nrow(inFileRdata)==0){Gdata2=1}
	 if(nrow(inFilegen2)==0 & nrow(inFileRdata)!=0){Gdata2=1}
	 if(nrow(inFilegen2)==0 & nrow(inFileRdata)==0){Gdata2=0}
	 
	 validate(
      need(Gdata2 != 0, "Please select second data file")	  
    )
	
	if(input$startAna=="StarChrom"){		
		if(SelFile=="Data"){
			
			withProgress(message = 'Getting...', value = 0,{
			incProgress(1/2, detail = "Wait, Please!")
			
			dirfile=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$datapath)    
			posit=data.frame(GenInfo$posit)
			#if(input$startAna=="StarBio"){
			#	posit=data.frame(GenInfo$posit)
			#}else{
			#	posit=vcfR::read.vcfR(as.character(dirfile))
			#	posit=data.frame(posit@fix)
			#}
			vcf1 <- vcfR::read.vcfR(as.character(inFilegen2$datapath))
			posit1=data.frame(vcf1@fix)
			filename1=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$name)
			filenamef=str_replace(filename1,".csv","")
			filename1=strsplit(filename1,"\\.")[[1]][1]	
			filename2=as.character(parseFilePaths(roots=getVolumes(), input$fileVCF)$name)
			filename2=strsplit(filename2,"\\.")[[1]][1]
			
			outFolder <- cambia_caracter(paste0("ChromMap_",filename1))
			setwd(str_replace(dirfile,filenamef,""))
			if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
			setwd("Output_BIO-R")
			if(!file.exists(outFolder)) dir.create(outFolder)
			setwd(outFolder)
			
			frames<-CompChrom(posit,posit1,filename1,filename2)
			chromosome_file=data.frame(frames[[1]])
			annotation_file=data.frame(frames[[2]])
			annotation_file2=data.frame(frames[[3]])	
			save(chromosome_file,annotation_file,annotation_file2,file="MapChrom.RData")
			incProgress(1, detail = "Finish")
			Sys.sleep(1)
			})
		}else{
			withProgress(message = 'Getting...', value = 0,{
			incProgress(1/2, detail = "Wait, Please!")
				load(as.character(inFileRdata$datapath))
			incProgress(1, detail = "Finish")
			Sys.sleep(1)
			})
		}
	}
	
	chromoMap(list(chromosome_file,chromosome_file),list(annotation_file,annotation_file2),n_win.factor = 5,export.options=T, fixed.window=F, plot.shift=T, remove.last.window=T, id="chrom2",ploidy=2,anno_col = c("green","red"))
	
  })

##################################################################################################################################################################  
  ### RUN BIOR- Diversity#########################################################################################################################
##################################################################################################################################################################
  DoforDiv<-reactive({
	if(input$startAna=="StarBio"){go1=1}else{go1=0}
	 validate(      
	  need(go1 != 0, "Please select the correspond type analysis (Biodiversity)")		  
    )
    outFolder<-"BioAnalysis"
    nall=2
    
	if(input$typedata=="vcfile"){
		datos<-as.data.frame(GenInfo$dfgen)
		headerdatos <- colnames(datos)
		#colnames(datos)=headerdatos
		newcolnames <- cambia_caracter(quita_espacio(as.matrix(colnames(datos))))
		colnames(datos) <- putg(newcolnames)
		groupfile=cbind(BioGID=colnames(datos)[-1],OriginalGID=headerdatos[-1])
	}else{
		if(input$typedata=="CUENTA"){
			datos<-GenInfo$dfgen
			headerdatos <- colnames(datos)
			newcolnames <- cambia_caracter(quita_espacio(as.matrix(headerdatos)))
			colnames(datos) <- putg(newcolnames)
			groupfile=cbind(BioGID=colnames(datos)[-1],OriginalGID=headerdatos[-1])
		}else{
			datos<-as.data.frame(GenInfo$dfgen)
			#inFilegen=parseFilePaths(roots=getVolumes(), input$filegen)
			headerdatos <- colnames(datos)
			#headerdatos <- read.table(as.character(inFilegen$datapath),nrows=1,header = FALSE, sep =',', stringsAsFactors = FALSE)
			#colnames(datos)=headerdatos
			newcolnames <- cambia_caracter(quita_espacio(as.matrix(headerdatos)))
			colnames(datos) <- putg(newcolnames)
			groupfile=cbind(BioGID=colnames(datos)[-1],OriginalGID=headerdatos[-1])
			#groupfile=cbind(BioGID=colnames(datos)[-1],OriginalGID=t(as.matrix(headerdatos[-1]))[,1])
		}
	}
	
    outFolder <- cambia_caracter(outFolder)
    distk<- cambia_caracter(input$distk)
    typedata<- cambia_caracter(input$typedata)
    missval=as.numeric(input$missval)
    mayorque=as.numeric(input$mayorque)
    menorque=as.numeric(input$menorque)
    ht1=as.numeric(input$ht1)
    ht2=as.numeric(input$ht2)
    ht3=as.numeric(input$ht3)
    
	if (typedata=="SNP"){
		datos=replace(datos,datos==ht1,99)
		datos=replace(datos,datos==ht2,0.5)
		datos=replace(datos,datos==ht3,999)
		datos=replace(datos,datos==99,1)   
		datos=replace(datos,datos==999,0)   
		datos <- data.frame(datos)
    }
	if (typedata=="vcfile"){
		datos=replace(datos,datos==0,99)
		datos=replace(datos,datos==1,0.5)
		datos=replace(datos,datos==2,999)
		datos=replace(datos,datos==99,1)   
		datos=replace(datos,datos==999,0)   
		datos <- data.frame(datos)
    }
    if (typedata=="FREQ"){datos <- data.frame(datos)}
	    
    req(input$filegen)
    dirfile=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$datapath)
    filename=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$name)
	filename1=strsplit(filename,"\\.")[[1]][1]	
	
    outFolder <- cambia_caracter(paste0("DiversityAnalysis_",filename1))
	setwd(str_replace(dirfile,filename,""))
    if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
    setwd("Output_BIO-R")
    if(!file.exists(outFolder)) dir.create(outFolder)
    setwd(outFolder)
    write.csv(groupfile,"ForGroups.csv",row.names=F, quote=F)
	system2('open',args="ForGroups.csv",wait=F)
	
	### Correr funciones ----------------------------------------
    if (typedata=="CUENTA"){
		#dtmp=as.data.frame(cbind(Allele=datos[,1],SNP=as.numeric(rep(c(NA,2)),dim(datos)[1]/2),datos[,-1]))
		Allele <- fread(as.character(dirfile), select = c(1))
		datos <- CounToFreqNI(datos,Allele,newcolnames)
	}
	if (typedata=="DistMat"){
		mrdMAT=as.matrix(datos[,-1])
		rownames(mrdMAT)=colnames(datos)[-1]
		colnames(mrdMAT)=colnames(datos)[-1]		
		##graphical representation, the MDS (multidimensional scaling) analysis
		mds=cmdscale(mrdMAT, k=3,eig=T)
		coord=mds$points
		perctCP12=c(round(mds$eig[1]/sum(mds$eig)*100,2),round(mds$eig[2]/sum(mds$eig)*100,2),round(mds$eig[3]/sum(mds$eig)*100,2))
		colnames(coord)=c("dim1","dim2","dim3")
		rm(mds)
		#########################################################################
		#########################################################################
		###For do dendograms and MDSgraph
		#########################################################################
		#########################################################################
		if (length(colnames(mrdMAT))<30){maxK=length(colnames(mrdMAT))-1} else {maxK=30}
		if (input$gapS==TRUE){
		#if (dim(mrdMAT)[1]<500){
		ver=fviz_nbclust(mrdMAT, hcut, nstart = 25, k.max=maxK,method = "gap_stat", nboot = 100)
		gapk1=ver[["data"]]$gap-ver[["data"]]$SE.sim
		test=cbind(ver[["data"]]$gap[-maxK],gapk1[-1])
		BestNc=min(which(test[,1]>=test[,2]))
		if (BestNc==1) {BestNc=2; print("Optimization fail (k=1)")}
		if (is.infinite(BestNc)==T) {BestNc=2; print("Optimization fail (k=InF)")}
		}else{
			BestNc=3
			#print("Optimization is not performed because there are many individuals and it takes a long time.")
		}
		clust=agnes(mrdMAT, method = "ward")
		coord2=cbind(gen=colnames(datos)[-1],coord)
		names(coord2)=c("Gen","Factor1","Factor2","Factor3")
		datos=NULL
		div=NULL		
		biodata=list(as.data.frame(div),coord2, getwd(), clust, datos, mrdMAT, perctCP12,BestNc)
	}else{
		biodata=Biodv(str_replace(filename,".csv",""),datos,nall,distk,mayorque,menorque,missval,typedata,ht1,ht2,ht3,input$gapS)		
	}
	updateTextInput(session,'tx','X Axis Label',value = paste0('Factor 1 (',biodata[[7]][1],'%)'))
	updateTextInput(session,'ty','Y Axis Label',value = paste0('Factor 2 (',biodata[[7]][2],'%)'))
	updateTextInput(session,'tz','Z Axis Label',value = paste0('Factor 3 (',biodata[[7]][3],'%)'))
	updateTextInput(session,'tx3D','X Axis Label',value = paste0('Factor 1 (',biodata[[7]][1],'%)'))
	updateTextInput(session,'ty3D','Y Axis Label',value = paste0('Factor 2 (',biodata[[7]][2],'%)'))
	updateTextInput(session,'tz3D','Z Axis Label',value = paste0('Factor 3 (',biodata[[7]][3],'%)'))
	updateTextInput(session,'nclust','No. Clusters',value = paste0(biodata[[8]]))
	
	return(biodata)	
  })

##################################################################################################################################################################  
  #Ver datos en tabla dinamica Summary Diversity
##################################################################################################################################################################
  output$seeDataDiver<-DT::renderDataTable({
	
	 if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	
	if(SelFile=="Data"){
		if(input$typedata=="DistMat"){
			seedatos=as.data.frame("Option no available for distance matrix input file")
		}else{	
			filename=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$name)
			file_name=str_replace(filename,".csv","")
			markfile=paste(DoforDiv()[[3]],"/MarkersRep_",file_name,".csv",sep="")
			genfile=paste(DoforDiv()[[3]],"/GenotypesRep_",file_name,".csv",sep="")
			if (file.exists(markfile) && file.exists(genfile)){
			shinyalert(title = "Important message", 
						text="Repeated sequence of markers and genotypes information were found!",
						type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
			)
			}
			if (file.exists(markfile) && !file.exists(genfile)){
			shinyalert(title = "Important message", 
						text="Repeated sequence of markers information were found!",
						type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
			)
			}
			if (!file.exists(markfile) && file.exists(genfile)){
			shinyalert(title = "Important message", 
						text="Repeated sequence of genotypes information were found!",
						type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
			)
			}
			
			seedatos<-as.data.frame(DoforDiv()[[1]])
		}
	}else{
		if(input$typedata=="DistMat"){
			seedatos=as.data.frame("Option no available for distance matrix input file")
		}else{
			UpRD=GenInfo$UploadRd
			seedatos<-UpRD$DivAna[[1]]
		}
	}
    datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
    
  })
  
  DoforPopStr<-reactive({
	 if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)	 
 	if(SelFile=="Data"){	
		datos<-as.data.frame(DoforDiv()[[5]])
		datos1<-as.data.frame(mdata1()[[1]])
		setwd(DoforDiv()[[3]])
		if(!file.exists("Output_MarkMonoGroups")) dir.create("Output_MarkMonoGroups")
		withProgress(message = 'Getting...', value = 0,{
		incProgress(1/2, detail = "Wait, Please!")
			if(typeof(datos1[,input$catv])!="double"){
				PopStr=gdiv(datos,datos1,input$catv,input$quitomono,as.data.frame(DoforDiv()[[6]]))
			}else{
				dt1=data.frame(rbind(t(c("Option no available","Option no available for continuous variables")),
                     t(c("for continuous variables","Option no available for continuous variables"))
				))
				PopStr=list(dt1,dt1)
			}
		incProgress(1, detail = "Finish")
		Sys.sleep(1)
		})
	}else{
		UpRD=GenInfo$UploadRd
		datos<-UpRD$DivAna[[5]]
		datos1<-UpRD$Aux[[1]]		
		setwd(as.character(UpRD$DivAna[[3]]))
		if(!file.exists("Output_MarkMonoGroups")) dir.create("Output_MarkMonoGroups")
		withProgress(message = 'Getting...', value = 0,{
		incProgress(1/2, detail = "Wait, Please!")
			if(typeof(datos1[,input$catv])!="double"){
				PopStr=gdiv(datos,datos1,input$catv,input$quitomono,as.data.frame(UpRD$DivAna[[6]]))
			}else{
				dt1=data.frame(rbind(t(c("Option no available","Option no available for continuous variables")),
                     t(c("for continuous variables","Option no available for continuous variables"))
				))
				PopStr=list(dt1,dt1)
			}
		incProgress(1, detail = "Finish")
		Sys.sleep(1)
		})
	}
	
	return(PopStr)
  })
  
##################################################################################################################################################################  
  #Ver datos en tabla dinamica Population structure
##################################################################################################################################################################
  output$seeDataGDiver<-DT::renderDataTable({
  #req(DoforPopStr())
	 if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	
	if(input$typedata=="DistMat"){
		seedatos=as.data.frame("Option no available for distance matrix input file")
	}else{
		seedatos=as.data.frame(DoforPopStr()[[1]])
	}	
	datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
  })
##################################################################################################################################################################  
  #Ver datos en tabla dinamica AMOVA
##################################################################################################################################################################
  output$seeDataGAmova<-DT::renderDataTable({
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	req(mdata1())
	if(SelFile=="Data"){	  
		if(input$typedata=="DistMat"){
			datos1<-as.data.frame(mdata1()[[1]])
			if(typeof(datos1[,input$catv])!="double"){
				pp=as.data.frame(as.character(datos1[,input$catv]))	
				rownames(pp)=datos1[,1]       
				groups=pp        
				agc.env=as.data.frame(as.numeric(as.factor(groups[,1])))
				names(agc.env)<-c("Pop")
				agc.env$Pop<-as.factor(agc.env$Pop)
				seedatos=forAMOVA(as.data.frame(DoforDiv()[[6]]),agc.env)  
				rownames(seedatos)=seedatos[,1]
				seedatos=seedatos[,-1]	
				write.csv(seedatos,file.path(DoforDiv()[[3]],paste0("AMOVA_",input$catv,".csv")))
			}else{
				seedatos=as.data.frame("Option no available for continuous variables")
			}
		}else{
			seedatos=as.data.frame(DoforPopStr()[[2]]) 
			rownames(seedatos)=seedatos[,1]
			seedatos=data.frame(seedatos[,-1])
		}
	}else{
		seedatos=as.data.frame(DoforPopStr()[[2]]) 
		rownames(seedatos)=seedatos[,1]
		seedatos=data.frame(seedatos[,-1])
	}
	datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
  })
##################################################################################################################################################################  
  #Transformacion de los datos, para el uso posterior en los graficos
  shinyFileChoose(input, 'fileenvbio', roots = getVolumes(),filetypes=c('', 'csv'))
##################################################################################################################################################################
  mdata1=reactive({
  if(input$startAna=="StarBio"){go1=1}else{go1=0}
	 validate(      
	  need(go1 != 0, "Please select the correspond type analysis (Biodiversity)")		  
    )
    #Cada que se actualice nclust
    pp= as.data.frame(cutree (as.hclust(DoforDiv()[[4]]), k = input$nclust))
    TFArx=as.phylo(as.hclust(DoforDiv()[[4]]))
    groups=as.data.frame(pp)
    coord2=as.data.frame(DoforDiv()[[2]])
    data1=as.data.frame(cbind(coord2,groups[,1]))
    names(data1)=c("Gen","Factor1","Factor2","Factor3","GroupClust")
    data1$Factor1=as.numeric(as.character(data1$Factor1))
    data1$Factor2=as.numeric(as.character(data1$Factor2))
    data1$Factor3=as.numeric(as.character(data1$Factor3))
    data1$GroupClust=as.factor(as.character(data1$GroupClust))
    
    #Cuando se agrega un archivo para grupos externos
    checkfile<-nrow(parseFilePaths(roots=getVolumes(), input$fileenvbio))
    if (checkfile!=0){
    inFileenvbio<-parseFilePaths(roots=getVolumes(), input$fileenvbio)
    dfenvbio <- read.csv(as.character(inFileenvbio$datapath),header = TRUE,sep = ",")
	dfenvbio[,1]<-putg(cambia_caracter(quita_espacio(as.character(dfenvbio[,1]))))
	indexCOV <- match(data1$Gen,as.character(dfenvbio[,1]))
    if(length(indexCOV)>0)	dfenvbio <- dfenvbio[indexCOV,]  
    #dfenvbio[,2] <- as.factor(dfenvbio[,2])	
    usenames=names(dfenvbio)[-1]
	dfenvbio<-dfenvbio[match(data1$Gen,dfenvbio[,1]),]
    for(i in 1:dim(dfenvbio)[2]){if (typeof(dfenvbio[,i])!="double"){dfenvbio[,i]=as.factor(as.character(dfenvbio[,i]))}}
    data1<-cbind(data1,dfenvbio[,-1])	
	names(data1)=c("Gen","Factor1","Factor2","Factor3","GroupClust",usenames)
    }
    #Nombre de las variables en la base de datos
    vars=names(data1)
    #Actualiza el selectinput de acuerdo a la base de datos cargada
    updateSelectInput(session,'xcol', 'X Variable',choices = vars,selected=vars[2])
    updateSelectInput(session,'ycol', 'Y Variable',choices = vars,selected=vars[3])
    updateSelectInput(session,'zcol', 'Z Variable',choices = vars,selected=vars[4])
    updateSelectInput(session,'catv', 'Group',choices = vars,selected=vars[5])
    updateSelectInput(session,'eti', 'Label',choices = vars,selected=vars[1])
    updateSelectInput(session,'xcol3D', 'X Variable',choices = vars,selected=vars[2])
    updateSelectInput(session,'ycol3D', 'Y Variable',choices = vars,selected=vars[3])
    updateSelectInput(session,'zcol3D', 'Z Variable',choices = vars,selected=vars[4])
    updateSelectInput(session,'catv3D', 'Group',choices = vars,selected=vars[5])
    updateSelectInput(session,'eti3D', 'Label',choices = vars,selected=vars[1])
    updateSelectInput(session,'catvdend', 'Group',choices = vars,selected=vars[5])
    
    result=list(data1,TFArx)
	write.csv(data1,file.path(DoforDiv()[[3]],"GroupsCluster.csv"),row.names=F)
    return(result)
  })
  #Actualiza la variable catv (grupos) en conjunto con la seleccion de colores
  #Cada que se elige un grupo se genera una cantidad de colores correspondiente 
  #al numero de elementos de cada grupo
  observeEvent(input$catv,{
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	set.seed(7)
    colores=colors()[-c(1,3:12,13:25,24,37:46,57:67,80,82,83,85:89,101:106,108:113,126:127,138,140:141,152:253,260:366,377:392,
                        394:447,449,478:489,492,513:534,536:546,557:561,579:583,589:609,620:629,418,436,646:651)]
    if(SelFile=="RData"){
		UpRD=GenInfo$UploadRd
		if(typeof(UpRD$Aux[[1]][,input$catv])!="double"){
			d=sample(colores,100)
			var=as.factor(UpRD$Aux[[1]][,input$catv])
			grupos=nlevels(var)
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}else{
		if(typeof(mdata1()[[1]][,input$catv])!="double"){
			d=sample(colores,100)
			var=as.factor(mdata1()[[1]][,input$catv])
			grupos=nlevels(var)
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}
    updateSelectInput(session,'color','Choose a color',choices=d,selected=d[1:grupos])
  })
  observeEvent(input$catv3D,{
    if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	set.seed(7)
    colores=colors()[-c(1,3:12,13:25,24,37:46,57:67,80,82,83,85:89,101:106,108:113,126:127,138,140:141,152:253,260:366,377:392,
                        394:447,449,478:489,492,513:534,536:546,557:561,579:583,589:609,620:629,418,436,646:651)]
    #d=sample(colores,100)
	if(SelFile=="RData"){
		UpRD=GenInfo$UploadRd
		if(typeof(UpRD$Aux[[1]][,input$catv3D])!="double"){
			d=sample(colores,100)
			var=as.factor(UpRD$Aux[[1]][,input$catv3D])
			grupos=nlevels(var)
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}else{
		if(typeof(mdata1()[[1]][,input$catv3D])!="double"){
			d=sample(colores,100)
			var=as.factor(mdata1()[[1]][,input$catv3D])
			grupos=nlevels(var)
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}
    updateSelectInput(session,'color3D','Choose a color', choices=d,selected=d[1:grupos])
  })
  observeEvent(input$catvdend,{
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
    set.seed(7)
    colores=colors()[-c(1,3:12,13:25,24,37:46,57:67,80,82,83,85:89,101:106,108:113,126:127,138,140:141,152:253,260:366,377:392,
                        394:447,449,478:489,492,513:534,536:546,557:561,579:583,589:609,620:629,418,436,646:651)]    
	#d=sample(colores,100)
	if(SelFile=="RData"){
		UpRD=GenInfo$UploadRd
		if(typeof(UpRD$Aux[[1]][,input$catvdend])!="double"){
			d=sample(colores,100)
			var=as.factor(UpRD$Aux[[1]][,input$catvdend])
			grupos=nlevels(var)
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}else{
		if(typeof(mdata1()[[1]][,input$catvdend])!="double"){
			d=sample(colores,100)
			var=as.factor(mdata1()[[1]][,input$catvdend])
			grupos=nlevels(var)
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}
    updateSelectInput(session,'colordend','Choose a color', choices=d,selected=d[1:grupos])
  })
##################################################################################################################################################################
  #grafico 2d
##################################################################################################################################################################
  output$try=renderPlotly({
	 if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	  
		if(!file.exists("Output_2DPlots")) dir.create("Output_2DPlots")    
		if(typeof(mdata1()[[1]][,input$catv])!="double"){
			p=plot_ly(data=mdata1()[[1]],x=mdata1()[[1]][,input$xcol],y=mdata1()[[1]][,input$ycol],color=mdata1()[[1]][,input$catv],
					type="scatter",mode="markers",colors = input$color,xaxis=F, yaxis=F,
					text=mdata1()[[1]][,input$eti],marker=list(size=input$size))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgp)%>%
			#titulo y etiquetas ejes
			layout(title=input$tp,titlefont=list(size=input$ts,color=input$pnc), xaxis = list(title = input$tx, titlefont=list(size=input$szl,color=input$ac)),
					yaxis = list(title = input$ty,titlefont=list(size=input$szl,color=input$ac)))
		}else{
			p=plot_ly(data=mdata1()[[1]],x=mdata1()[[1]][,input$xcol],y=mdata1()[[1]][,input$ycol],color=mdata1()[[1]][,input$catv],
					type="scatter",mode="markers",xaxis=F, yaxis=F,colors=c("blue","cyan","green","orange","red"),
					text=mdata1()[[1]][,input$eti],marker=list(size=input$size))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgp)%>%
			#titulo y etiquetas ejes
			layout(title=input$tp,titlefont=list(size=input$ts,color=input$pnc), xaxis = list(title = input$tx, titlefont=list(size=input$szl,color=input$ac)),
					yaxis = list(title = input$ty,titlefont=list(size=input$szl,color=input$ac)))		
		}
		#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 2d
		#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
		withr::with_dir(file.path(DoforDiv()[[3]],"Output_2DPlots"),saveWidget(p,paste0('MDS2d_',input$catv,'.html'), selfcontained = F))
		p
	}else{
		UpRD=GenInfo$UploadRd
		if(!file.exists("Output_2DPlots")) dir.create("Output_2DPlots")    
		if(typeof(UpRD$Aux[[1]][,input$catv])!="double"){
			p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcol],y=UpRD$Aux[[1]][,input$ycol],color=UpRD$Aux[[1]][,input$catv],
					type="scatter",mode="markers",colors = input$color,xaxis=F, yaxis=F,
					text=UpRD$Aux[[1]][,input$eti],marker=list(size=input$size))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgp)%>%
			#titulo y etiquetas ejes
			layout(title=input$tp,titlefont=list(size=input$ts,color=input$pnc), xaxis = list(title = input$tx, titlefont=list(size=input$szl,color=input$ac)),
					yaxis = list(title = input$ty,titlefont=list(size=input$szl,color=input$ac)))
		}else{
			p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcol],y=UpRD$Aux[[1]][,input$ycol],color=UpRD$Aux[[1]][,input$catv],
					type="scatter",mode="markers",xaxis=F, yaxis=F,colors=c("blue","cyan","green","orange","red"),
					text=UpRD$Aux[[1]][,input$eti],marker=list(size=input$size))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgp)%>%
			#titulo y etiquetas ejes
			layout(title=input$tp,titlefont=list(size=input$ts,color=input$pnc), xaxis = list(title = input$tx, titlefont=list(size=input$szl,color=input$ac)),
					yaxis = list(title = input$ty,titlefont=list(size=input$szl,color=input$ac)))
		}
		#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 2d
		#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
		withr::with_dir(file.path(UpRD$DivAna[[3]],"Output_2DPlots"),saveWidget(p,paste0('MDS2d_',input$catv,'.html'), selfcontained = F))
		p
	}
  })
  
  #Cuadro de texto que muestra la direccion en la que se guardo el grafico.
  output$default1=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	  
		paste('You can find results and edited plots files in:',as.character(DoforDiv()[[3]]))
	}else{
		UpRD=GenInfo$UploadRd
		paste('You can find results and edited plots files in:',as.character(UpRD$DivAna[[3]]))
	}
  })
##################################################################################################################################################################  
  #grafico 3d
##################################################################################################################################################################
  output$try3d=renderPlotly({
    if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	 
		if(!file.exists("Output_3DPlots")) dir.create("Output_3DPlots")
		if(typeof(mdata1()[[1]][,input$catv3D])!="double"){	
			p=plot_ly(data=mdata1()[[1]],x=mdata1()[[1]][,input$xcol3D],y=mdata1()[[1]][,input$ycol3D],z=mdata1()[[1]][,input$zcol3D],color=mdata1()[[1]][,input$catv3D],
					type = 'scatter3d' ,mode="markers",colors = input$color3D,
					text=mdata1()[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
			#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
					paper_bgcolor=input$bkgp3D,
					title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		}else{
			p=plot_ly(data=mdata1()[[1]],x=mdata1()[[1]][,input$xcol3D],y=mdata1()[[1]][,input$ycol3D],z=mdata1()[[1]][,input$zcol3D],color=mdata1()[[1]][,input$catv3D],
				type = 'scatter3d' ,mode="markers",colors=c("blue","cyan","green","orange","red"),
				text=mdata1()[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
		#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
							yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
							zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
				paper_bgcolor=input$bkgp3D,
				title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		
		}
		#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 3d
		#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
		withr::with_dir(file.path(DoforDiv()[[3]],"Output_3DPlots"),saveWidget(p,paste0('MDS3d_',input$catv3D,'.html'), selfcontained = F))
		p
	}else{
		UpRD=GenInfo$UploadRd
		if(!file.exists("Output_3DPlots")) dir.create("Output_3DPlots")
		if(typeof(UpRD$Aux[[1]][,input$catv3D])!="double"){		
			p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcol3D],y=UpRD$Aux[[1]][,input$ycol3D],z=UpRD$Aux[[1]][,input$zcol3D],color=UpRD$Aux[[1]][,input$catv3D],
					type = 'scatter3d' ,mode="markers",colors = input$color3D,
					text=UpRD$Aux[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
			#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
					paper_bgcolor=input$bkgp3D,
					title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		}else{
			p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcol3D],y=UpRD$Aux[[1]][,input$ycol3D],z=UpRD$Aux[[1]][,input$zcol3D],color=UpRD$Aux[[1]][,input$catv3D],
					type = 'scatter3d' ,mode="markers",colors=c("blue","cyan","green","orange","red"),
					text=UpRD$Aux[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
			#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
					paper_bgcolor=input$bkgp3D,
					title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		}
		#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 3d
		#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
		withr::with_dir(file.path(UpRD$DivAna[[3]],"Output_3DPlots"),saveWidget(p,paste0('MDS3d_',input$catv3D,'.html'), selfcontained = F))
		p
	}
  })
  
  #Cuadro de texto que muestra la direccion en la que se guardo el grafico.
  output$default3d=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	  
	  paste('You can find results and edited plots files in:',as.character(DoforDiv()[[3]]))
  }else{
	  UpRD=GenInfo$UploadRd
	  paste('You can find results and edited plots files in:',as.character(UpRD$DivAna[[3]]))
  }
 })
##################################################################################################################################################################
#distance matrix heatmap
##################################################################################################################################################################
  output$heat=renderPlotly({  
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	req(mdata1())
	if(SelFile=="Data"){	  
		use=as.data.frame(mdata1()[[1]])
		group=input$catv
		use=use[order(use[,group]),]
		useorder=match(use$Gen,rownames(DoforDiv()[[6]]))
		distplot=plot_ly(x=rownames(DoforDiv()[[6]])[useorder],y=rownames(DoforDiv()[[6]])[useorder],z = DoforDiv()[[6]][useorder,useorder], colorscale=input$colorheat,type = "heatmap")%>%
		layout(xaxis = list(showticklabels = F), yaxis = list(showticklabels = F))
		withr::with_dir(file.path(DoforDiv()[[3]]),saveWidget(distplot,'DistancesPlotEdited.html', selfcontained = F))
		distplot
	}else{
		UpRD=GenInfo$UploadRd
		use=as.data.frame(UpRD$Aux[[1]])
		group=input$catv
		use=use[order(use[,group]),]
		useorder=match(use$Gen,rownames(UpRD$DivAna[[6]]))
		distplot=plot_ly(x=rownames(UpRD$DivAna[[6]])[useorder],y=rownames(UpRD$DivAna[[6]])[useorder],z = UpRD$DivAna[[6]][useorder,useorder], colorscale=input$colorheat,type = "heatmap")%>%
		layout(xaxis = list(showticklabels = F), yaxis = list(showticklabels = F))
		withr::with_dir(file.path(UpRD$DivAna[[3]]),saveWidget(distplot,'DistancesPlotEdited.html', selfcontained = F))
		distplot
	}
  })
  
  output$defaultheat=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	  
		paste('You can find results and edited plots files in:',as.character(DoforDiv()[[3]]))
	}else{
		UpRD=GenInfo$UploadRd
		paste('You can find results and edited plots files in:',as.character(UpRD$DivAna[[3]]))
	}
  })
##################################################################################################################################################################
#dendogram plot
##################################################################################################################################################################  
  output$dend=renderPlot({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	
		DivAna=DoforDiv()
		Aux=mdata1()
		save(DivAna,Aux,file="DivAna.RData")
		if(!file.exists("Output_Dendograms")) dir.create("Output_Dendograms")
		data=as.data.frame(mdata1()[[1]])
		info<- data[,c("Gen",input$catvdend)]
		info<- cbind(ID=info$Gen,info)
		names(info)=c("ID","Gen","Group")
		tree=mdata1()[[2]]
		if(typeof(info$Group)!="double"){
			p=ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
			scale_color_manual(values=input$colordend)+
			geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, hjust=0.5)+ 
			theme(legend.position=input$poslen)
		}else{
			p=ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
			scale_color_gradientn(colours=c("blue","cyan","green","orange","red")) +
			geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, hjust=0.5)+ 
			theme(legend.position=input$poslen)
		}
		ggsave(paste0(DoforDiv()[[3]],'\\Output_Dendograms\\DendogramPlot_',input$catvdend,'.pdf'),p)
		p
	}else{
		UpRD=GenInfo$UploadRd
		if(!file.exists("Output_Dendograms")) dir.create("Output_Dendograms")
		data=as.data.frame(UpRD$Aux[[1]])
		info<- data[,c("Gen",input$catvdend)]
		info<- cbind(ID=info$Gen,info)
		names(info)=c("ID","Gen","Group")
		tree=UpRD$Aux[[2]]
		if(typeof(info$Group)!="double"){
			p=ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
			scale_color_manual(values=input$colordend)+
			geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, hjust=0.5)+ 
			theme(legend.position=input$poslen)
		}else{
			p=ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
			scale_color_gradientn(colours=c("blue","cyan","green","orange","red")) +
			geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, hjust=0.5)+ 
			theme(legend.position=input$poslen)
		}
		ggsave(paste0(UpRD$DivAna[[3]],'\\Output_Dendograms\\DendogramPlot_',input$catvdend,'.pdf'),p)
		p
	}
  })
  
  
  output$defaultdend=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 validate(
      need(Gdata != 0, "Please upload data")
	)
	if(SelFile=="Data"){	
		paste('You can find results and edited plots files in:',as.character(DoforDiv()[[3]]))
	}else{
		UpRD=GenInfo$UploadRd
		paste('You can find results and edited plots files in:',as.character(UpRD$DivAna[[3]]))
	}  
  })
  
  ########################################################################################################################################################################################    
  #Core-BIO-R
  ########################################################################################################################################################################################    
  shinyFileChoose(input, 'filedistbio', roots = getVolumes(),filetypes=c('', 'csv'))
  shinyFileChoose(input, 'filephendatbio', roots = getVolumes(),filetypes=c('', 'csv'))
  
  observe({
  
    if("phendat"%in%input$datause){	  
      shinyjs::show(id="filephendatbio")
	  shinyjs::show(id="gdEN")
      shinyjs::show(id="gdAN")
      shinyjs::show(id="gdEE")
    }else{
      shinyjs::hide(id="filephendatbio")
	  shinyjs::hide(id="gdEN")
      shinyjs::hide(id="gdAN")
      shinyjs::hide(id="gdEE")
    }
	
	if("distdat"%in%input$datause){	 
      shinyjs::show(id="filedistbio")
	  shinyjs::show(id="pcdEN")
      shinyjs::show(id="pcdAN")
      shinyjs::show(id="pcdEE")
    }else{
      shinyjs::hide(id="filedistbio")
	  shinyjs::hide(id="pcdEN")
      shinyjs::hide(id="pcdAN")
      shinyjs::hide(id="pcdEE")
    }
	
	if("gendat"%in%input$datause){		
	  shinyjs::show(id="mrdEN")
      shinyjs::show(id="csedEN")
      shinyjs::show(id="mrdAN")
      shinyjs::show(id="csedAN")
      shinyjs::show(id="mrdEE")
      shinyjs::show(id="csedEE")
      shinyjs::show(id="SH")
      shinyjs::show(id="HE")
      shinyjs::show(id="CV")
	}else{
	  shinyjs::hide(id="mrdEN")
      shinyjs::hide(id="csedEN")
      shinyjs::hide(id="mrdAN")
      shinyjs::hide(id="csedAN")
      shinyjs::hide(id="mrdEE")
      shinyjs::hide(id="csedEE")
      shinyjs::hide(id="SH")
      shinyjs::hide(id="HE")
      shinyjs::hide(id="CV")
	}
	
	if("gendat"%in%input$datause && "distdat"%in%input$datause){
	  shinyjs::hide(id="filedistbio")
	  shinyjs::hide(id="mrdEN")
      shinyjs::hide(id="csedEN")
      shinyjs::hide(id="mrdAN")
      shinyjs::hide(id="csedAN")
      shinyjs::hide(id="mrdEE")
      shinyjs::hide(id="csedEE")
      shinyjs::hide(id="SH")
      shinyjs::hide(id="HE")
      shinyjs::hide(id="CV")
	  shinyjs::hide(id="pcdEN")
      shinyjs::hide(id="pcdAN")
      shinyjs::hide(id="pcdEE")
	}
	
	if("phendat"%in%input$datause && "distdat"%in%input$datause){		
	  shinyjs::hide(id="filedistbio")
	  shinyjs::hide(id="filephendatbio")
	  shinyjs::hide(id="gdEN")
      shinyjs::hide(id="gdAN")
      shinyjs::hide(id="gdEE")
	  shinyjs::hide(id="pcdEN")
      shinyjs::hide(id="pcdAN")
      shinyjs::hide(id="pcdEE")
	}
	
	if("gendat"%in%input$datause && "distdat"%in%input$datause && "phendat"%in%input$datause){
	  shinyjs::show(id="filedistbio")
	  shinyjs::show(id="filephendatbio")
	  shinyjs::show(id="gdEN")
      shinyjs::show(id="gdAN")
      shinyjs::show(id="gdEE")
	  shinyjs::show(id="pcdEN")
      shinyjs::show(id="pcdAN")
      shinyjs::show(id="pcdEE")
	  shinyjs::show(id="mrdEN")
      shinyjs::show(id="csedEN")
      shinyjs::show(id="mrdAN")
      shinyjs::show(id="csedAN")
      shinyjs::show(id="mrdEE")
      shinyjs::show(id="csedEE")
      shinyjs::show(id="SH")
      shinyjs::show(id="HE")
      shinyjs::show(id="CV")
	}
	
  })
  
  
  datacore=reactive({
    if ("phendat"%in%input$datause){
      typedata<- cambia_caracter(input$typedata)
	  if(!"gendat"%in%input$datause && "distdat"%in%input$datause){checadist="exist"}else{checadist="none"}
	  
	  dirfilePhen<-parseFilePaths(roots=getVolumes(), input$filephendatbio)
      validate(
		need(nrow(dirfilePhen)!=0, "Please select data"),
		need(checadist!="exist","Combination of data NO available")
      )
    }else{      
      dirfilePhen<-"none"      
    }
  
    if ("gendat"%in%input$datause){
      datosgen=as.data.frame(GenInfo$dfgen)
	  #inFilegen=parseFilePaths(roots=getVolumes(), input$filegen)
	  #headerdatos <- read.table(as.character(inFilegen$datapath),nrows=1,header = FALSE, sep =',', stringsAsFactors = FALSE)
      #colnames(datosgen)=headerdatos    
      dirfileGen<-"exist"
	  typedata<- cambia_caracter(input$typedata)	  
	  validate(
        need(nrow(datosgen) != 0, "Please select data"),
		need(typedata!="DistMat","Please upload genetic data")
      )      
      ht1=as.numeric(input$ht1)
      ht2=as.numeric(input$ht2)
      ht3=as.numeric(input$ht3)
      newcolnames <- cambia_caracter(quita_espacio(as.matrix(colnames(datosgen))))
      colnames(datosgen) <- putg(newcolnames)
	  if (typedata=="vcfile"){
		datosgen=replace(datosgen,datosgen==0,99)
		datosgen=replace(datosgen,datosgen==1,0.5)
		datosgen=replace(datosgen,datosgen==2,999)
		datosgen=replace(datosgen,datosgen==99,1)   
		datosgen=replace(datosgen,datosgen==999,0)   
		datosgen <- as.data.frame(datosgen)
	  }
      if (typedata=="SNP"){
        datosgen=replace(datosgen,datosgen==ht1,99)
		datosgen=replace(datosgen,datosgen==ht2,0.5)
		datosgen=replace(datosgen,datosgen==ht3,999)
		datosgen=replace(datosgen,datosgen==99,1)   
		datosgen=replace(datosgen,datosgen==999,0)   
		datosgen <- as.data.frame(datosgen)
      }
      if (typedata=="FREQ"){datosgen <-as.data.frame(datosgen)}
	  
	  if (typedata=="CUENTA"){
		dtmp=as.data.frame(cbind(Allele=datosgen[,1],SNP=as.numeric(rep(c(NA,2)),dim(datosgen)[1]/2),datosgen[,-1]))
		datosgen <- as.data.frame(CounToFreqNI(data_set=dtmp))
	  }
	  Marker=unlist(lapply(seq(1:dim(datosgen)[1]),function(y) rbind(y,y)))
      Allele=rep(c("",".1"),dim(datosgen)[1])
	  MAlle=paste("M",Marker,Allele,sep="")
      mcomp=apply(datosgen[,-1],2,function(y) rbind(y,1-y))
	  NGen=colnames(mcomp)
	  mcomp=data.frame(t(mcomp))
	  datos=as.data.frame(cbind(NGen,mcomp))
	  colnames(datos)=c("NAME",MAlle)
	  rownames(datos)=seq(1:length(NGen))
	  rm(datosgen)
	  gc()
    }else{
      datos=as.data.frame(matrix(0,3,3))
      dirfileGen<-"none"
      
      } 
    
    if ("distdat"%in%input$datause){      
	  typedata<- cambia_caracter(input$typedata)	  
      dirfileDist<-parseFilePaths(roots=getVolumes(), input$filedistbio)
      validate(
        need(nrow(dirfileDist) != 0, "Please select data")
      )
    }else{
      dirfileDist<-"none" 
    }
    
    ver=sum(as.numeric(c(input$mrdEN,input$csedEN,input$gdEN,input$pcdEN,
              input$mrdEE,input$csedEE,input$gdEE,input$pcdEE,
              input$mrdAN,input$csedAN,input$gdAN,input$pcdAN,
              input$SH,input$HE,input$CV)))
    validate(
      need(ver == 1, "The sum of the weights must be one")
    )
    
    if("phendat"%in%input$datause){
      dirfile=as.character(parseFilePaths(roots=getVolumes(), input$filephendatbio)$datapath)
      filename=as.character(parseFilePaths(roots=getVolumes(), input$filephendatbio)$name)
    }else{
      if("gendat"%in%input$datause){
        dirfile=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$datapath)
        filename=as.character(parseFilePaths(roots=getVolumes(), input$filegen)$name)
      }else{
        dirfile=as.character(parseFilePaths(roots=getVolumes(), input$filedistbio)$datapath)
        filename=as.character(parseFilePaths(roots=getVolumes(), input$filedistbio)$name)
      }
    }
	outFolder <- cambia_caracter(paste("CoreSubset_",str_replace(filename,".csv",""),sep=""))
    setwd(str_replace(dirfile,filename,""))
    if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
    setwd("Output_BIO-R")
    if(!file.exists(outFolder)) dir.create(outFolder)
    setwd(outFolder)
    fin=getwd()
	
	### Correr funciones ----------------------------------------
	withProgress(message = 'Getting...', value = 0,{
	incProgress(1/2, detail = "Wait, Please!")
	
    funchunt(datos,dirfileGen,dirfilePhen,dirfileDist,input$score,input$mrdEN,input$csedEN,input$gdEN,
             input$pcdEN,input$mrdAN,input$csedAN,input$gdAN,input$pcdAN,input$mrdEE,input$csedEE,input$gdEE,input$pcdEE,
             input$SH,input$HE,input$CV)
    incProgress(1, detail = "Finish")
	Sys.sleep(1)
	})
	return(fin)
  })
  
  output$defaultcore=renderText({
    HTML(paste0("<font color=\"#FF0000\"><b> ","SUCCESSFUL ANALYSIS!!", "</b></font>","<br> You can find results files in:","<font color=\"#FF0000\"><b>",datacore(),"</b></font>"))
    })
##################################################################################################################################################################
###################################################################################################################################################################################################
#Start GWAS Code
##################################################################################################################################################################
##################################################################################################################################################################
shinyFileChoose(input, 'filephen', roots = getVolumes(),filetypes=c('csv'))
shinyFileChoose(input, 'fileplants', roots = getVolumes(),filetypes=c('gz'))

myDataPhen<-reactive({
	inFilephen=parseFilePaths(roots=getVolumes(), input$filephen)
	if(nrow(inFilephen)==0){PhenData=0}else{PhenData=1}	 
	validate(
		need(PhenData != 0, "Please select phen data")
    )
	dropsPheno<-fread(as.character(inFilephen$datapath),header = TRUE,sep = ",",na=c("NA",".","-",""))
	names(dropsPheno)[1]="genotype"
	
	datos=data.frame(GenInfo$dfgen)	
	c1=length(which(colnames(datos)[-1]%in%dropsPheno$genotype==F))
	validate(
		need(c1 == 0, "No match gen and pheno info, please check")
	)
	
	varsP=names(dropsPheno)[-1]
	updateSelectInput(session,'ManTrait', 'Trait',choices = varsP,selected=varsP[1])
	updateSelectInput(session,'QQTrait', 'Trait',choices = varsP,selected=varsP[1])
	
	return(list(dropsPheno))
})

myDataRegGen<-reactive({
	inFileplants=parseFilePaths(roots=getVolumes(), input$fileplants)
	if(nrow(inFileplants)==0){PlantsData=0}else{PlantsData=1}
	validate(
		need(PlantsData != 0, "Please select genome ref data")
    )
	print("Please wait...")
	dirdata=str_replace(parseFilePaths(roots=getVolumes(), input$filegen)$datapath,parseFilePaths(roots=getVolumes(), input$filegen)$name,"")
	load(paste0(dirdata,"Output_BIO-R\\ChromMap_",cambia_caracter(strsplit(parseFilePaths(roots=getVolumes(), input$filegen)$name,"[.]","")[[1]][1]),"\\MapChrom.RData"))	
	wheat<-import(as.character(inFileplants$datapath))
	wheat=as.data.frame(wheat)	
	check=sort(as.character(unique(wheat[,1])[1:dim(chromosome_file)[1]]))
	checknum=as.numeric(levels(as.factor(check)))
	wheat[,1]=as.character(wheat[,1])
	for (k in 1:dim(chromosome_file)[1]){
		wheat[which(wheat[,1]==check[k]),1]=checknum[k]
	}
	wheat=cbind(wheat[,c(1,12,7,2,3,8,5,9)],paste0("gene_id \"",wheat[,10],"\"; transcript_id \"",wheat[,13],"\"; exon_number \"",wheat[,16],"\";"))
	names(wheat)=c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
	wheat[,1]=as.character(wheat[,1])
	wheat[,2]=as.character(wheat[,2])
	wheat[,3]=as.character(wheat[,3])
	wheat[,6]=as.character(wheat[,6])
	wheat[,7]=as.character(wheat[,7])
	wheat[,8]=as.character(wheat[,8])
	wheat[,9]=as.character(wheat[,9])
	
	
	return(list(wheat))
})

DoforGWAS<-reactive({
		if(input$startAna=="StarChrom" & input$typedata=="vcfile"){
		req(myDataPhen())
		dirdata=str_replace(parseFilePaths(roots=getVolumes(), input$filegen)$datapath,parseFilePaths(roots=getVolumes(), input$filegen)$name,"")
		load(paste0(dirdata,"Output_BIO-R\\ChromMap_",cambia_caracter(strsplit(parseFilePaths(roots=getVolumes(), input$filegen)$name,"[.]","")[[1]][1]),"\\MapChrom.RData"))
		
		datos<-GenInfo$hapmap
		datos<-datos[,-c(2:11)]
		newcolnames <- cambia_caracter(quita_espacio(as.matrix(colnames(datos))))
		colnames(datos) <- putg(newcolnames)
		datos<-datos[which(datos[,1]%in%annotation_file$ID==T),]
		rownames(datos)<-datos[,1]
		datos<-t(datos[,-1])
		
		posit<-GenInfo$posit
		posit<-posit[,c(3,1,2)]
		posit<-posit[which(posit[,1]%in%annotation_file$ID==T),]
		names(posit)<-c("SNP.names","chr","pos")
		posit$chr=as.numeric(posit$chr)
		posit$pos=as.numeric(posit$pos)
		rownames(posit)=posit$SNP.names
		
		dropsPheno<-data.frame(myDataPhen()[[1]])
		newcolnames <- cambia_caracter(quita_espacio(as.matrix(dropsPheno$genotype)))
		dropsPheno$genotype<- putg(newcolnames)	 	
		
		gDataDrops <- createGData(geno = datos,map = posit,pheno = dropsPheno)
		gDataDrops <- codeMarkers(gData = gDataDrops,
                           nMissGeno = as.numeric(input$nmissgeno), 
                           nMiss = as.numeric(input$nmissind), 
                           MAF= as.numeric(input$nmf),
                           impute = TRUE, 
                           imputeType = "random",
                           verbose=TRUE)
		Ttraits=as.character(input$gwastraits)
		Ttraits=unlist(strsplit(Ttraits,","))
		GWASDrops <- runSingleTraitGwas(gData = gDataDrops, traits = Ttraits)
		
		outFolderPh <- cambia_caracter(paste0("GWASAnalysis_",strsplit(parseFilePaths(roots=getVolumes(), input$filegen)$name,"[.]","")[[1]][1]))
		setwd(dirdata)
		if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
		setwd("Output_BIO-R")
		if(!file.exists(outFolderPh)) dir.create(outFolderPh)
		setwd(outFolderPh)
		
		save(gDataDrops,GWASDrops,file="DataGWAS.RData")
		for(nt in 1:length(Ttraits)){			
			ass=split(GWASDrops[["GWAResult"]][["dropsPheno"]],GWASDrops[["GWAResult"]][["dropsPheno"]][["trait"]])
			ass=ass[Ttraits[nt]]
			write.csv(ass,paste0("MarkersEffect_",Ttraits[nt],".csv"))		
		}
	return(list(GWASDrops,savein=as.character(getwd())))
	}
})
	

DoforGene<-reactive({
	if(input$startAna=="StarChrom" & input$typedata=="vcfile"){
	req(myDataRegGen())
	wheat<-as.data.frame(myDataRegGen()[[1]])
	#chrdim<-as.data.frame(myDataRegGen()[[2]])
	GWASDrops<-DoforGWAS()[[1]]
	
	ass=split(GWASDrops[["GWAResult"]][["dropsPheno"]],GWASDrops[["GWAResult"]][["dropsPheno"]][["trait"]])
	ass=ass[[as.character(input$RegGenPlotraits)]]
	ass=ass[,c(2,3,4,6)]
	names(ass)=c("Marker","Locus","Site","p")
	
	#vcf <- vcfR::read.vcfR(parseFilePaths(roots=getVolumes(), input$filegen)$datapath)
	#htm=vcfR2hapmap(vcf)
	#htm<-data.frame(htm[-1,])
	htm<-GenInfo$hapmap
	htm$chrom=as.numeric(htm$chrom)
	htm$pos=as.numeric(htm$pos)
	
	savein<-DoforGWAS()[[2]]
	setwd(savein)
	save(GWASDrops,ass,wheat,htm,file="DataGWASRegGen.RData")
	
	return(list(ass,wheat,htm,savein=getwd()))
	}
})

##################################################################################################################################################################
#For see regional associated gene plot
##################################################################################################################################################################
  output$defaultRegGenPlotChr=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	
	if(SelFile=="Data"){	  
		paste('You can find results and edited plots files in:',as.character(DoforGene()[[4]]))
	}else{
		dirdata=str_replace(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath,parseFilePaths(roots=getVolumes(), input$fileRData)$name,"")
		paste('You can find results and edited plots files in:',as.character(dirdata))
	}
  })
  output$RegGenPlotChrPlot=renderPlot({
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	if(SelFile=="Data"){
		ass=DoforGene()[[1]]
		wheat=DoforGene()[[2]]
		htm=DoforGene()[[3]]
		#chrdim=DoforGene()[[5]]
		savein=DoforGene()[[4]]
	}else{
		load(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath)
	}
	
	wheatmp=wheat[which(wheat[,1]==as.character(input$RegGenPlotChr)),]
	
	IntReg=IntRegionalPlot(chr=as.numeric(input$RegGenPlotChr),left=as.numeric(input$RegGenPlotleft),
	right=as.numeric(input$RegGenPlotrigth),threshold=as.numeric(input$RegGenPlotChrThr),gtf=wheat,association=ass,
	hapmap=htm,label_gene_name=TRUE,hapmap_ld=htm,leadsnp_size=3)
	listgen=IntReg[["plot_env"]][["gene_list"]]
	
	setwd(savein)
	if(dim(listgen)[1]!=0){	write.csv(listgen,paste0("Genelist_",input$RegGenPlotraits,"_chr",input$RegGenPlotChr,".csv"))}
	pdf(paste0("RegionalPlot_",input$RegGenPlotraits,"_chr",input$RegGenPlotChr,".pdf"), width = 6.89, height = 5.41)
		print(IntReg)
	dev.off()
	print(IntReg)
  })
  
##################################################################################################################################################################
#For see manhattan plot
##################################################################################################################################################################
  output$defaultMan=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	if(SelFile=="Data"){	  
		paste('You can find results and edited plots files in:',as.character(DoforGWAS()[[2]]))
	}else{
		dirdata=str_replace(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath,parseFilePaths(roots=getVolumes(), input$fileRData)$name,"")
		paste('You can find results and edited plots files in:',as.character(dirdata))
	}
  })
  output$ManPlot=renderPlot({
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	if(SelFile=="Data"){
		GWASDrops=DoforGWAS()[[1]]
		savein=DoforGWAS()[[2]]
	}else{
		load(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath)
	}
	
	if(input$Manchr=="All"){ManchrV=1:length(unique(GWASDrops[["GWAResult"]][["dropsPheno"]][["chr"]]))}else{ManchrV=eval(parse(text=input$Manchr))}
	if(input$MancolPalette=="Blue"){MancolorPaletteV=rep(c("blue","darkblue"),length(ManchrV))[1:length(ManchrV)]}
	if(input$MancolPalette=="Gray"){MancolorPaletteV=rep(c("black","cornsilk4"),length(ManchrV))[1:length(ManchrV)]}
	if(input$MancolPalette=="Green"){MancolorPaletteV=rep(c("darkolivegreen1","darkolivegreen4"),length(ManchrV))[1:length(ManchrV)]}
	if(input$MancolPalette=="Colors"){MancolorPaletteV=colorRampPalette(brewer.pal(11,"Spectral"))(length(ManchrV))}
	
	plot(GWASDrops, plotType = "manhattan", trait = input$ManTrait, yThr = as.numeric(input$ManyThr) ,chr= ManchrV ,colPalette= MancolorPaletteV)
	
	setwd(savein)
	#jpeg(paste0("ManhathanPlot_",input$ManTrait,"_",input$ManyThr,".jpeg"), width = 800, height = 600)
	pdf(paste0("ManhathanPlot_",input$ManTrait,"_",input$ManyThr,".pdf"), width = 7, height = 6)
		plot(GWASDrops, plotType = "manhattan", trait = input$ManTrait, yThr = as.numeric(input$ManyThr) ,chr= ManchrV ,colPalette= MancolorPaletteV)
	dev.off()
  })
  
##################################################################################################################################################################
#For see QTL plot
##################################################################################################################################################################
  output$defaultQTL=renderText({
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	if(SelFile=="Data"){	  
		paste('You can find results and edited plots files in:',as.character(DoforGWAS()[[2]]))
	}else{
		dirdata=str_replace(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath,parseFilePaths(roots=getVolumes(), input$fileRData)$name,"")
		paste('You can find results and edited plots files in:',as.character(dirdata))
	}
  })
  output$QTLPlot=renderPlot({
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	if(SelFile=="Data"){
		GWASDrops=DoforGWAS()[[1]]
		savein=as.character(DoforGWAS()[[2]])
	}else{
		load(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath)
	}
	
	if(input$QTLchr=="All"){QTLchrV=1:length(unique(GWASDrops[["GWAResult"]][["dropsPheno"]][["chr"]]))}else{QTLchrV=eval(parse(text=input$QTLchr))}
	
	plot(GWASDrops, plotType = "qtl", yThr = as.numeric(input$QTLyThr) ,chr= QTLchrV )
	
	setwd(savein)
	#jpeg(paste0("QTLPlot_",input$QTLyThr,".jpeg"), width = 800, height = 600)
	pdf(paste0("QTLPlot_",input$QTLyThr,".pdf"), width = 7, height = 6.3)
		plot(GWASDrops, plotType = "qtl", yThr = as.numeric(input$QTLyThr) ,chr= QTLchrV )
	dev.off()
  })
  
##################################################################################################################################################################
#For see QQ plot
##################################################################################################################################################################
  output$defaultQQ=renderText({
	#req(myDataGen())
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	req(DoforGWAS())
	
	if(SelFile=="Data"){	  
		paste('You can find results and edited plots files in:',as.character(DoforGWAS()[[2]]))
	}else{
		dirdata=str_replace(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath,parseFilePaths(roots=getVolumes(), input$fileRData)$name,"")
		paste('You can find results and edited plots files in:',as.character(dirdata))
	}
  })
  output$QQPlot=renderPlot({
	#req(myDataGen())
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(input$startAna=="StarChrom" & input$typedata=="vcfile"){test=1}else{test=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(test != 0, "No option available, please select Chromosom option")
	)
	req(DoforGWAS())	
	
	if(SelFile=="Data"){
		GWASDrops=DoforGWAS()[[1]]
		savein=as.character(DoforGWAS()[[2]])
	}else{
		load(parseFilePaths(roots=getVolumes(), input$fileRData)$datapath)
	}
	
	if(input$QTLchr=="All"){QTLchrV=1:length(unique(GWASDrops[["GWAResult"]][["dropsPheno"]][["chr"]]))}else{QTLchrV=eval(parse(text=input$QTLchr))}
	#If the lambda value is greater than 1, then this may be evidence for some systematic bias that needs to be corrected in your analysis.
	print("Lambda for GWAS:")
	print(GWASDrops[["GWASInfo"]][["inflationFactor"]])
	
	plot(GWASDrops, plotType = "qq", trait=input$QQTrait)
	
	setwd(savein)
	#jpeg(paste0("QQPlot_",input$QQTrait,".jpeg"), width = 800, height = 600)
	pdf(paste0("QQPlot_",input$QQTrait,".pdf"), width = 6, height = 6)
		plot(GWASDrops, plotType = "qq", trait=input$QQTrait)
	dev.off()
  })



##################################################################################################################################################################
}



