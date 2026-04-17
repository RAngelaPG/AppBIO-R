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
  #shinyFileChoose(input, 'filegen', roots = getVolumes(),filetypes=c('', 'csv','vcf','gz'))
  #shinyFileChoose(input, 'fileRdata', roots = getVolumes(),filetypes=c('', 'RData'))
  
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
	
	if(input$startAna=="CompChrom"){
		shinyjs::show(id="fileVCF")
	}
	if(input$startAna=="StarChrom"){
		shinyjs::hide(id="typedata")
		shinyjs::hide(id="missval")
		shinyjs::hide(id="missvalG")
		shinyjs::hide(id="mayorque")
		shinyjs::hide(id="menorque")	
		shinyjs::hide(id="fileVCF")
		shinyjs::hide(id="ht1")
		shinyjs::hide(id="ht2")
		shinyjs::hide(id="ht3")
		updateRadioButtons(session, "typedata", selected = "vcfile")		
	}	
	if(input$startAna=="pavs"){
		shinyjs::hide(id="typedata")
		shinyjs::hide(id="missval")
		shinyjs::hide(id="missvalG")
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
		shinyjs::show(id="missvalG")
		shinyjs::show(id="mayorque")
		shinyjs::show(id="menorque")
		shinyjs::hide(id="fileVCF")		
	}
	if(input$gapS){shinyjs::runjs('document.getElementById("nclust").style.pointerEvents = "none";')}else{shinyjs::runjs('document.getElementById("nclust").style.pointerEvents = "auto";')}
	
  })
  
  
	GenInfo <- reactiveValues(dfgen=NULL,posit=NULL,UploadRd=NULL,hapmap=NULL,pavs=NULL,summaryChrom=NULL,chromosome_file=NULL,annotation_file=NULL, annotation_file2=NULL)
	observeEvent(input$readgenofile, {  
		UploadRd=NULL

	inFilegen=input$filegen	
	inFileRdata=input$fileRdata	
	 if(!is.null(inFilegen) & is.null(inFileRdata)){Gdata=1;SelFile="Data"}
	 if(is.null(inFilegen) & !is.null(inFileRdata)){Gdata=1;SelFile="RData"}
	 if(is.null(inFilegen) & is.null(inFileRdata)){Gdata=0}
	 if(!is.null(inFilegen) & !is.null(inFileRdata)){NAdata=0}else{NAdata=1} 
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
				pavs = fread(as.character(inFilegen$datapath),header = TRUE,sep = ",",na=c("NA",".","-",""))
				posit=NULL
				hapmap=NULL
				dfgen=NULL
				if(length(intersect(file_attr_PAVS, pavs[4,])) == 15){
					gename=as.matrix(pavs[4,16:dim(pavs)[2]])
					idname=as.character(as.matrix(pavs[-c(1:4),1]))
					pavs=as.data.frame(pavs[-c(1:4),-c(1:15)])
					pavs=apply(pavs,2,as.numeric)
					colnames(pavs)=gename
					rownames(pavs)=idname
				}else{
					gename=names(pavs)[-1]
					idname=as.character(as.matrix(pavs[,1]))
					pavs=as.data.frame(pavs[,-1])
					pavs=apply(pavs,2,as.numeric)
					colnames(pavs)=gename
					rownames(pavs)=idname
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
				pavs=NULL				
			}else{
			if(input$typedata=="CUENTA"){
				dfgen = read.csv.ffdf(file=as.character(inFilegen$datapath),VERBOSE=T)
				posit=NULL
				hapmap=NULL
				pavs=NULL
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
				pavs=NULL
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
			pavs=NULL
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
			#updateSelectInput(session,'catv3D', 'Group',choices = vars,selected=vars[5])
			updateSelectInput(session,'eti3D', 'Label',choices = vars,selected=vars[1])
			#updateSelectInput(session,'catvdend', 'Group',choices = vars,selected=vars[5])    
			updateSelectInput(session,'xcolCH', 'X Variable',choices = vars,selected=vars[2])
			updateSelectInput(session,'ycolCH', 'Y Variable',choices = vars,selected=vars[3])
			updateSelectInput(session,'zcolCH', 'Z Variable',choices = vars,selected=vars[4])
			updateSelectInput(session,'catvCH', 'Group',choices = vars,selected=vars[5])
			updateSelectInput(session,'etiCH', 'Label',choices = vars,selected=vars[1])			
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
		shinyjs::hide(id="methgap")
		shinyjs::hide(id="mixture")
		shinyjs::hide(id="nclust")
		shinyjs::hide(id="typedata")
		shinyjs::hide(id="distk")
	}else{
		shinyjs::show(id="fileenvbio")
		shinyjs::show(id="quitomono")
		shinyjs::show(id="gapS")
		shinyjs::show(id="methgap")
		shinyjs::show(id="mixture")		
		shinyjs::show(id="nclust")
		shinyjs::show(id="distk")
	}
	
    GenInfo$dfgen<-dfgen
	GenInfo$posit<-posit
	GenInfo$UploadRd<-UploadRd
	GenInfo$hapmap<-hapmap
	GenInfo$pavs<-pavs
	
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
	hasData   <- !is.null(GenInfo$pavs)
	hasRData  <- !is.null(GenInfo$UploadRd)
	bothData  <- hasData && hasRData
	anyData   <- hasData || hasRData
	correctTab <- input$startAna == "pavs"

	validate(
		need(anyData, "Please upload data"),
		need(!bothData, "Both data selected, please close and try again"),
		need(correctTab, "Go to the correct analysis tab")
	)
 
  # seleccionar fuente
  if (hasData) {
    use <- GenInfo$pavs
  } else {
    e <- new.env()
    load(input$fileRdata$datapath, envir = e)
    use <- e$GenInfo$pavs
  }  
  
  pp<-plotly::plot_ly(x = colnames(use), y = rownames(use), z = use,  type = "heatmap",  colorscale = "Jet" ) %>%
    layout( xaxis = list(showticklabels = TRUE), yaxis = list(showticklabels = FALSE))
  pp	
		
  })
 
 output$download_pavs <- downloadHandler(
  filename = function() {"PAVSPlot.html"},
  content = function(file) {
	req(GenInfo$pavs)
    use=GenInfo$pavs
    pp<-plotly::plot_ly(x=colnames(use),y=rownames(use),z = use, colorscale="Jet",type = "heatmap")%>%
			layout(xaxis = list(showticklabels = T), yaxis = list(showticklabels = F))	    
    saveWidget(pp, file, selfcontained = T)
  }
)

##################################################################################################################################################################
  #Ver grafico de cromosomas
##################################################################################################################################################################
  tablaChrom<-reactive({	
	inFilegen=input$filegen
	inFileRdata=input$fileRdata	
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
			
		if(SelFile=="Data"){
			#dirfile=as.character(input$filegen$datapath)    		
			posit=data.frame(GenInfo$posit)				
			tmark=dim(posit)[1]
			posit$POS=as.numeric(posit$POS)
			posit$CHROM=as.character(posit$CHROM)
			posit$ID=as.character(posit$ID)
			tabamv1=data.frame(table(as.factor(posit$CHROM)))				
			colnames(tabamv1)=c("Chromosome","NMarks")
			dimt1=dim(tabamv1)[1]
			
			chromi=as.character(tabamv1[which(tabamv1$NMarks>=1),1])				
			chromi1=chromi[grep("chr",chromi)]
			if(length(chromi1)==0){chromi=chromi}else{chromi=chromi1}
			if(length(which(posit$CHROM%in%chromi==F))!=0){posit$CHROM[which(posit$CHROM%in%chromi==F)]=NA}		
			if(length(which(is.na(posit[,1])))!=0){posit=posit[-which(is.na(posit[,1])),]}
			tabamv1=data.frame(table(as.factor(posit$CHROM)))		
			colnames(tabamv1)=c("V1","V2")
			#colnames(tabamv1)=c("Chromosome","NMarks")
			if(dimt1>dim(tabamv1)[1]) {
				shinyalert(title = "Important message", 
					text="Extra chromosome information was \n found in the blast and was deleted",closeOnEsc = FALSE, 
					type = "warning", showCancelButton=FALSE, showConfirmButton=TRUE, confirmButtonCol = "green"
				)
					#print("Extra chromosome information was found in the blast and was deleted") 
			}			
					
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
			names(chrdim)=c("V1","V2")			
			tabamv1=rbind(rbind(c("Total markers: ",tmark),c("Markers in blast: ",dim(posit)[1]),c("Markers in blast for each chromosome: "," "),c("Chromosome","NMarks")),tabamv1)						
			#write.csv(chrdim,"ChromDim.csv")
		    GenInfo$chromosome_file=chromosome_file
			GenInfo$annotation_file=annotation_file
			#save(chromosome_file,annotation_file,file="MapChrom.RData")
			GenInfo$summaryChrom=rbind(tabamv1,c("minPos","MaxPos"),chrdim)
		}else{
			load(as.character(inFileRdata$datapath))			
		}			
			
	}	
	#chromoMap(list(chromosome_file),list(annotation_file),n_win.factor = 5, export.options=T, fixed.window=F, remove.last.window=T,plot.shift=T, id="chrom1")	
  })
  
  output$seeChromPlot<-renderChromoMap({	
	inFilegen=input$filegen
	inFileRdata=input$fileRdata	
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
	
    tablaChrom()	
	chromoMap(list(GenInfo$chromosome_file),list(GenInfo$annotation_file),n_win.factor = 5, export.options=T, fixed.window=F, remove.last.window=T,plot.shift=T, id="chrom1")
	
	})
	
  output$tablaChrom <- DT::renderDataTable({
   req(GenInfo$summaryChrom)
	withProgress(message = 'Getting...', value = 0,{
		incProgress(1/2, detail = "Wait, Please!")	
		#tablaChrom()
		if(input$startAna=="StarChrom"){summd<-as.data.frame(GenInfo$summaryChrom)}else{summd<-NULL}
		incProgress(1, detail = "Finish")
		Sys.sleep(1)
	})
	nombre <- paste0("SummaryChromoMap_", Sys.Date())
	datatable(summd, selection="multiple", escape=FALSE, extensions = 'Buttons',
              options = list(dom  = 'Bfrtip',buttons=list(list(extend='csv',filename=nombre),list(extend='excel',filename=nombre)),pageLength = 7,width="100%", scrollX = TRUE))		  
	
  })
  #shinyFileChoose(input, 'fileVCF', roots = getVolumes(),filetypes=c('vcf','gz'))

##################################################################################################################################################################
  #Ver grafico comparacion de blast cromosomas
##################################################################################################################################################################
  BtablaChrom<-reactive({	

	useData  <-  !is.null(GenInfo$dfgen)
	useRData <- !is.null(GenInfo$UploadRd)
	SelFile <- NULL	
	inFileVCF   <- input$fileVCF
	
	if (useData) { SelFile <- "Data"} else if (useRData) {  SelFile <- "RData"}		
	
		if(SelFile=="Data"){
			posit <- data.frame(GenInfo$posit)
			vcf1 <- vcfR::read.vcfR(inFileVCF$datapath)
			posit1 <- data.frame(vcf1@fix)

			filename1 <- tools::file_path_sans_ext(inFileVCF$name)
			filename2 <- tools::file_path_sans_ext(inFileVCF$name)

			frames <- CompChrom(posit, posit1, filename1, filename2)			
			GenInfo$chromosome_file  <- data.frame(frames[[1]])
			GenInfo$annotation_file  <- data.frame(frames[[2]])
			GenInfo$annotation_file2 <- data.frame(frames[[3]])
			GenInfo$summaryChrom     <- data.frame(frames[[4]])
			
		} else {
			e <- new.env()
			load(input$fileRdata$datapath, envir = e)
		}	
	
	 
  })
  
output$BChromPlot<-renderChromoMap({		
	inFilegen=input$filegen
	inFileRdata=input$fileRdata	
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 if(!is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){NAdata=0}else{NAdata=1}
	 if(input$typedata=="vcfile"){ dataT1=1 }else {dataT1=0}	 
	 if(input$startAna=="CompChrom"){go1=1}else{go1=0}
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(NAdata != 0, "Both data selected, please close and try again"),
	  need(dataT1 != 0, "Option only available for VCF files"),		  	  	  
	  need(go1 != 0, "go to the correspond tab of type analysis")		  
    )
	
		BtablaChrom()	        					
		chromoMap(list(GenInfo$chromosome_file,GenInfo$chromosome_file),list(GenInfo$annotation_file,GenInfo$annotation_file2),n_win.factor = 5,export.options=T, fixed.window=F, remove.last.window=T, plot.shift=T, id="chrom2",ploidy=2,anno_col = c("green","red"))
	
	})
	
  output$BtablaChrom <-DT::renderDataTable({
	req(GenInfo$summaryChrom)	
	withProgress(message = 'Getting...', value = 0,{
		incProgress(1/2, detail = "Wait, Please!")
		#BtablaChrom()
		if(input$startAna=="CompChrom"){summd2<-as.data.frame(GenInfo$summaryChrom)}else{summd2<-NULL}
		incProgress(1, detail = "Finish")
		Sys.sleep(1)
	})
	nombre <- paste0("SummaryChromoMapCompare_", Sys.Date())
	datatable(summd2, selection="multiple", escape=FALSE, extensions = 'Buttons',
              options = list(dom  = 'Bfrtip',buttons=list(list(extend='csv',filename=nombre),list(extend='excel',filename=nombre)),pageLength = 7,width="100%", scrollX = TRUE))		  
  })  
  
##################################################################################################################################################################  
  ### RUN BIOR- Diversity#########################################################################################################################
##################################################################################################################################################################

data_gen <- reactive({
  req(GenInfo$dfgen)	
		datos<-as.data.frame(GenInfo$dfgen)
		req(ncol(datos) > 1)
		headerdatos <- colnames(datos)		
		newcolnames <- cambia_caracter(quita_espacio(as.matrix(colnames(datos))))
		colnames(datos) <- putg(newcolnames)
		groupfile=cbind(BioGID=colnames(datos)[-1],OriginalGID=headerdatos[-1])	
	groupfile
})

output$downloadData <- downloadHandler(
  filename = function() {"ForGroups.csv"},
  content = function(file) {
    req(data_gen())
    write.csv(data_gen(), file, row.names = FALSE)
  }
)

 BiodivInfo <- reactiveValues(res1=NULL,res2=NULL)
	observeEvent(input$calcopt, {
  shinybusy::show_modal_spinner(spin = "fading-circle", text = "Running diversity analysis..." )
  #withProgress(message = "Running analysis...", value = 0, {

    #incProgress(0.3, detail = "Processing diversity...")
    BiodivInfo$res1 <- isolate(DoforDiv())

    #incProgress(0.6, detail = "Processing metadata...")
    BiodivInfo$res2 <- isolate(mdata1())

    #incProgress(1, detail = "Done")

  #}) 
	shinybusy::remove_modal_spinner()
})

	
  DoforDiv<-reactive({
	if(input$startAna=="StarBio"){go1=1}else{go1=0}
	 validate(      
	  need(go1 != 0, "Please select the correspond type analysis (Biodiversity)")		  
    )
    outFolder<-"BioAnalysis"
    nall=2
    	
		datos<-as.data.frame(GenInfo$dfgen)
		headerdatos <- colnames(datos)	
		newcolnames <- cambia_caracter(quita_espacio(as.matrix(colnames(datos))))
		colnames(datos) <- putg(newcolnames)
		#groupfile=cbind(BioGID=colnames(datos)[-1],OriginalGID=headerdatos[-1])
		
    #outFolder <- cambia_caracter(outFolder)
    distk<- cambia_caracter(input$distk)
    typedata<- cambia_caracter(input$typedata)
    missval=as.numeric(input$missval)
	missvalG=as.numeric(input$missvalG)
    mayorque=as.numeric(input$mayorque)
    menorque=as.numeric(input$menorque)
	mixture=input$mixture
	gap=input$gapS
	methodgap=input$methgap
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
    #dirfile=as.character(input$filegen$datapath)
    filename=as.character(input$filegen$name)
	#filename1=strsplit(filename,"\\.")[[1]][1]		
    #outFolder <- cambia_caracter(paste0("DiversityAnalysis_",filename1))
	#setwd(str_replace(dirfile,filename,""))
    #if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
    #setwd("Output_BIO-R")
    #if(!file.exists(outFolder)) dir.create(outFolder)
    #setwd(outFolder)
    #write.csv(groupfile,"ForGroups.csv",row.names=F, quote=F)
	#system2('open',args="ForGroups.csv",wait=F)
	#downloadHandler(filename = "ForGroups.csv", content = function(file) { write.csv(groupfile, file, row.names = FALSE)})	
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
		
		if(input$gapS==TRUE){	
		print("Do optimization...")
		if(input$methgap=="gap"){			
			ver=fviz_nbclust(mrdMAT, hcut, k.max=20,method = "gap_stat", nboot = 100)
			gapk1=ver[["data"]]$gap-ver[["data"]]$SE.sim
			test=cbind(ver[["data"]]$gap[-(dim(mrdMAT)[1]-2)],gapk1[-1])
			BestNc=which(test[,1]>=test[,2])  
			if(length(BestNc)>1){ if(min(BestNc)==1){BestNc=BestNc[2]}else{BestNc=min(BestNc)} }			
			if (BestNc==1) {BestNc=2; print("Optimization fail (k=1)")}
			if (is.infinite(BestNc)==T){BestNc=2; print("Optimization fail (k=InF)")}
		}else{
			ver <- fviz_nbclust(mrdMAT,hcut, method = "silhouette",k.max = 20)
			BestNc <- ver$data$clusters[which.max(ver$data$y)]
			if(length(BestNc)>1){ if(min(BestNc)==1){BestNc=BestNc[2]}else{BestNc=min(BestNc)} }
			if (BestNc==1) {BestNc=2; print("Optimization fail (k=1)")}
			if (is.infinite(BestNc)==T){BestNc=2; print("Optimization fail (k=InF)")}
		}		
		}else{BestNc=3}
	
		clust=agnes(mrdMAT, method = "ward")
		coord2=cbind(gen=colnames(datos)[-1],coord)
		names(coord2)=c("Gen","Factor1","Factor2","Factor3")
		datos=NULL
		div=NULL		
		biodata=list(as.data.frame(div),coord2, getwd(), clust, datos, mrdMAT, perctCP12,BestNc,qmatrix=NULL,exadiv=NULL, exadivg=NULL, tablaOpt=NULL)
	}else{
		biodata=Biodv(str_replace(filename,".csv",""),datos,nall,distk,mayorque,menorque,missval,typedata,ht1,ht2,ht3,missvalG,mixture,gap,methodgap,as.numeric(input$nclust))		
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
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	
	if(SelFile=="Data"){
		if(input$typedata=="DistMat"){
			seedatos=as.data.frame("Option no available for distance matrix input file")
		}else{				
			seedatos<-as.data.frame(BiodivInfo$res1[[1]])
			seedatos[,2]=round(as.numeric(seedatos[,2]),4)
			names(seedatos)=c("source","value")
		}
	}else{
		if(input$typedata=="DistMat"){
			seedatos=as.data.frame("Option no available for distance matrix input file")
		}else{
			UpRD=GenInfo$UploadRd
			seedatos<-UpRD$DivAna[[1]]
			seedatos[,2]=round(as.numeric(seedatos[,2]),4)
			names(seedatos)=c("source","value")
		}
	}
    datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
    
  })
  
  DoforPopStr<-reactive({
	 if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	) 
 	if(SelFile=="Data"){	
		datos<-as.data.frame(BiodivInfo$res1[[5]])
		datos1<-as.data.frame(BiodivInfo$res2[[1]])
		#setwd(BiodivInfo$res1[[3]])
		#if(!file.exists("Output_MarkMonoGroups")) dir.create("Output_MarkMonoGroups")
		withProgress(message = 'Getting...', value = 0,{
		incProgress(1/2, detail = "Wait, Please!")
			if(typeof(datos1[,input$catv])!="double"){
				PopStr=gdiv(datos,datos1,input$catv,input$quitomono,as.data.frame(BiodivInfo$res1[[6]]))
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
		#setwd(as.character(UpRD$DivAna[[3]]))
		#if(!file.exists("Output_MarkMonoGroups")) dir.create("Output_MarkMonoGroups")
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
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	
	if(input$typedata=="DistMat"){
		seedatos=as.data.frame("Option no available for distance matrix input file")
	}else{
		seedatos=as.data.frame(DoforPopStr()[[1]])
		seedatos[,3] <- ifelse(grepl("^\\d+(\\.\\d+)?$", seedatos[,3]),round(as.numeric(seedatos[,3]), 4), seedatos[,3])
		names(seedatos)=c("source","group","value","usefulMark")
	}	
	datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
  })
##################################################################################################################################################################  
  #Ver datos en tabla dinamica AMOVA
##################################################################################################################################################################
  amovaTable<-reactive({  
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	#req(mdata1())
	if(SelFile=="Data"){	  
		if(input$typedata=="DistMat"){
			datos1<-as.data.frame(BiodivInfo$res2[[1]])
			if(typeof(datos1[,input$catv])!="double"){
				pp=as.data.frame(as.character(datos1[,input$catv]))	
				rownames(pp)=datos1[,1]       
				groups=pp        
				agc.env=as.data.frame(as.numeric(as.factor(groups[,1])))
				names(agc.env)<-c("Pop")
				agc.env$Pop<-as.factor(agc.env$Pop)
				seedatos=forAMOVA(as.data.frame(BiodivInfo$res1[[6]]),agc.env)  
				rownames(seedatos)=seedatos[,1]
				seedatos=seedatos[,-1]	
				seedatos[,2:8]<-apply(seedatos[,2:8], 2,function(x) {round(as.numeric(x),4)})
				#write.csv(seedatos,file.path(BiodivInfo$res1[[3]],paste0("AMOVA_",input$catv,".csv")))
			}else{
				seedatos=as.data.frame("Option no available for continuous variables")
			}
		}else{
			seedatos=as.data.frame(DoforPopStr()[[2]]) 
			rownames(seedatos)=seedatos[,1]
			seedatos=data.frame(seedatos[,-1])
			seedatos[,2:8]<-apply(seedatos[,2:8], 2,function(x) {round(as.numeric(x),4)})
		}
	}else{
		seedatos=as.data.frame(DoforPopStr()[[2]]) 
		rownames(seedatos)=seedatos[,1]
		seedatos=data.frame(seedatos[,-1])
		seedatos[,2:8]<-apply(seedatos[,2:8], 2,function(x) {round(as.numeric(x),4)})
	}
	seedatos
  })
  
  output$seeDataGAmova<-DT::renderDataTable({
    seedatos=as.data.frame(amovaTable())
    datatable(seedatos, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
  })
##################################################################################################################################################################  
  #Transformacion de los datos, para el uso posterior en los graficos
  #shinyFileChoose(input, 'fileenvbio', roots = getVolumes(),filetypes=c('', 'csv'))
##################################################################################################################################################################
colores_18<-c("intense red" = "#E41A1C",
  "electric blue" = "#0057FF",
  "bright green" = "#00C853",
  "vivid yellow" = "#FFD600",
  "strong orange" = "#FF3D00",
  "intense purple" = "#AA00FF",
  "brown" = "#A65628",
  "pink" = "#F781BF",
  "cyan" = "#17BECF",
  "olive" = "#BCBD22",
  "dark blue" = "#08306B",
  "dark wine" = "#67000D",
  "forest green" = "#00441B",
  "deep purple" = "#3F007D",
  "dark brown" = "#4A2C2A",
  "bright sky blue" = "#00B0FF",
  "coral pink" = "#FF6F91",
  "dark lime green" = "#7CB342"
  )
  
  mdata1=reactive({
  if(input$startAna=="StarBio"){go1=1}else{go1=0}
	 validate(      
	  need(go1 != 0, "Please select the correspond type analysis (Biodiversity)")		  
    )	
	
    #Cada que se actualice nclust
	pp= as.data.frame(cutree (as.hclust(BiodivInfo$res1[[4]]), k = as.numeric(BiodivInfo$res1[[8]])))	
    TFArx=as.phylo(as.hclust(BiodivInfo$res1[[4]]))
    groups=as.data.frame(pp)
    coord2=as.data.frame(BiodivInfo$res1[[2]])
    data1=as.data.frame(cbind(coord2,groups[,1]))
    names(data1)=c("Gen","Factor1","Factor2","Factor3","GroupClust")
    data1$Factor1=as.numeric(as.character(data1$Factor1))
    data1$Factor2=as.numeric(as.character(data1$Factor2))
    data1$Factor3=as.numeric(as.character(data1$Factor3))
    data1$GroupClust=as.factor(as.character(data1$GroupClust))
    
    #Cuando se agrega un archivo para grupos externos
    #checkfile<-is.null(input$fileenvbio)
    if (!is.null(input$fileenvbio)){
    inFileenvbio<-input$fileenvbio
    dfenvbio <- read.csv(as.character(inFileenvbio$datapath),header = TRUE,sep = ",")
	if(!all(grepl("^g", dfenvbio[,1]))){
		dfenvbio[,1]<-putg(cambia_caracter(quita_espacio(as.character(dfenvbio[,1]))))
	}
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
	vgg=vars[5:length(vars)]
	vxyz=vars[2:4]
    #Actualiza el selectinput de acuerdo a la base de datos cargada
    updateSelectInput(session,'xcol', 'X Variable',choices = vxyz,selected=vxyz[1])
    updateSelectInput(session,'ycol', 'Y Variable',choices = vxyz,selected=vxyz[2])
    updateSelectInput(session,'zcol', 'Z Variable',choices = vxyz,selected=vxyz[3])
    updateSelectInput(session,'catv', 'Group',choices = vgg,selected=vgg[1])
    updateSelectInput(session,'eti', 'Label',choices = vars,selected=vars[1])
    updateSelectInput(session,'xcol3D', 'X Variable',choices = vxyz,selected=vxyz[1])
    updateSelectInput(session,'ycol3D', 'Y Variable',choices = vxyz,selected=vxyz[2])
    updateSelectInput(session,'zcol3D', 'Z Variable',choices = vxyz,selected=vxyz[3])
    updateSelectInput(session,'eti3D', 'Label',choices = vars,selected=vars[1])    
    
    result=list(data1,TFArx)
	#write.csv(data1,file.path(BiodivInfo$res1[[3]],"GroupsCluster.csv"),row.names=F)
    return(result)
  })
  #Actualiza la variable catv (grupos) en conjunto con la seleccion de colores
  #Cada que se elige un grupo se genera una cantidad de colores correspondiente 
  #al numero de elementos de cada grupo
  observeEvent(input$catv,{
	# Asegura que exista la variable seleccionada
    req(input$catv)
  
	if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	#hasResult <- !is.null(BiodivInfo$res1)
	  # Validaciones
	req(Gdata != 0)
	req(!is.null(BiodivInfo$res1))
  
	set.seed(7)
    colores=colors()[-c(1,3:12,13:25,24,37:46,57:67,80,82,83,85:89,101:106,108:113,126:127,138,140:141,152:253,260:366,377:392,
                        394:447,449,478:489,492,513:534,536:546,557:561,579:583,589:609,620:629,418,436,646:651)]
    if(SelFile=="RData"){
		UpRD=GenInfo$UploadRd
		if(typeof(UpRD$Aux[[1]][,input$catv])!="double"){			
			var=as.factor(UpRD$Aux[[1]][,input$catv])
			grupos=nlevels(var)
			if(grupos>18){d=sample(colores,100)}else{d=colores_18}
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}else{
		if(typeof(BiodivInfo$res2[[1]][,input$catv])!="double"){			
			var=as.factor(BiodivInfo$res2[[1]][,input$catv])
			grupos=nlevels(var)
			if(grupos>18){d=sample(colores,100)}else{d=colores_18}
		}else{
			d=c("Jet","Jet","Jet","Jet","Jet","Jet")
			grupos=2
		}
	}
    updateSelectInput(session,'color','Choose a color',choices=d,selected=d[1:grupos])
	updateSelectInput(session,'color3D','Choose a color', choices=d,selected=d[1:grupos])
	updateSelectInput(session,'colordend','Choose a color', choices=d,selected=d[1:grupos])	
  })
    
##################################################################################################################################################################
  #grafico 2d
##################################################################################################################################################################
  mds2Plot<-reactive({
	 if(!is.null(GenInfo$dfgen) && is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) && !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) && is.null(GenInfo$UploadRd)){Gdata=0}
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	if(SelFile=="Data"){	  
		#if(!file.exists("Output_2DPlots")) dir.create("Output_2DPlots")    
		if(typeof(BiodivInfo$res2[[1]][,input$catv])!="double"){
			p=plot_ly(data=BiodivInfo$res2[[1]],x=BiodivInfo$res2[[1]][,input$xcol],y=BiodivInfo$res2[[1]][,input$ycol],color=BiodivInfo$res2[[1]][,input$catv],
					type="scatter",mode="markers",colors = input$color,xaxis=F, yaxis=F,
					text=BiodivInfo$res2[[1]][,input$eti],marker=list(size=input$size))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgp)%>%
			#titulo y etiquetas ejes
			layout(title=input$tp,titlefont=list(size=input$ts,color=input$pnc), xaxis = list(title = input$tx, titlefont=list(size=input$szl,color=input$ac)),
					yaxis = list(title = input$ty,titlefont=list(size=input$szl,color=input$ac)))
		}else{
			p=plot_ly(data=BiodivInfo$res2[[1]],x=BiodivInfo$res2[[1]][,input$xcol],y=BiodivInfo$res2[[1]][,input$ycol],color=BiodivInfo$res2[[1]][,input$catv],
					type="scatter",mode="markers",xaxis=F, yaxis=F,colors=c("blue","cyan","green","orange","red"),
					text=BiodivInfo$res2[[1]][,input$eti],marker=list(size=input$size))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgp)%>%
			#titulo y etiquetas ejes
			layout(title=input$tp,titlefont=list(size=input$ts,color=input$pnc), xaxis = list(title = input$tx, titlefont=list(size=input$szl,color=input$ac)),
					yaxis = list(title = input$ty,titlefont=list(size=input$szl,color=input$ac)))		
		}
		#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 2d
		#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
		#withr::with_dir(file.path(BiodivInfo$res1[[3]],"Output_2DPlots"),saveWidget(p,paste0('MDS2d_',input$catv,'.html'), selfcontained = F))
		p
	}else{
		UpRD=GenInfo$UploadRd
		#if(!file.exists("Output_2DPlots")) dir.create("Output_2DPlots")    
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
		#withr::with_dir(file.path(UpRD$DivAna[[3]],"Output_2DPlots"),saveWidget(p,paste0('MDS2d_',input$catv,'.html'), selfcontained = F))
		p
	}
  })
  
 output$try=renderPlotly({
	mds2Plot()
 })
##################################################################################################################################################################  
  #grafico 3d
##################################################################################################################################################################
  mds3Plot<-reactive({
    if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	if(SelFile=="Data"){	 
		#if(!file.exists("Output_3DPlots")) dir.create("Output_3DPlots")
		if(typeof(BiodivInfo$res2[[1]][,input$catv])!="double"){	
			p=plot_ly(data=BiodivInfo$res2[[1]],x=BiodivInfo$res2[[1]][,input$xcol3D],y=BiodivInfo$res2[[1]][,input$ycol3D],z=BiodivInfo$res2[[1]][,input$zcol3D],color=BiodivInfo$res2[[1]][,input$catv],
					type = 'scatter3d' ,mode="markers",colors = input$color3D,
					text=BiodivInfo$res2[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
			#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
					paper_bgcolor=input$bkgp3D,
					title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		}else{
			p=plot_ly(data=BiodivInfo$res2[[1]],x=BiodivInfo$res2[[1]][,input$xcol3D],y=BiodivInfo$res2[[1]][,input$ycol3D],z=BiodivInfo$res2[[1]][,input$zcol3D],color=BiodivInfo$res2[[1]][,input$catv],
				type = 'scatter3d' ,mode="markers",colors=c("blue","cyan","green","orange","red"),
				text=BiodivInfo$res2[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
		#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
							yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
							zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
				paper_bgcolor=input$bkgp3D,
				title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		
		}
		#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 3d
		#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
		#withr::with_dir(file.path(BiodivInfo$res1[[3]],"Output_3DPlots"),saveWidget(p,paste0('MDS3d_',input$catv3D,'.html'), selfcontained = F))
		p
	}else{
		UpRD=GenInfo$UploadRd
		#if(!file.exists("Output_3DPlots")) dir.create("Output_3DPlots")
		if(typeof(UpRD$Aux[[1]][,input$catv])!="double"){		
			p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcol3D],y=UpRD$Aux[[1]][,input$ycol3D],z=UpRD$Aux[[1]][,input$zcol3D],color=UpRD$Aux[[1]][,input$catv],
					type = 'scatter3d' ,mode="markers",colors = input$color3D,
					text=UpRD$Aux[[1]][,input$eti3D],marker=list(size=input$size3D))%>%
			#Nombre de los ejes del grafico
			layout(scene=list(xaxis = list(title = input$tx3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								yaxis = list(title = input$ty3D,titlefont=list(size=input$szl3D,color=input$ac3D)),
								zaxis = list(title = input$tz3D,titlefont=list(size=input$szl3D,color=input$ac3D))),
					paper_bgcolor=input$bkgp3D,
					title=input$tp3D,titlefont=list(size=input$ts3D,color=input$pnc3D))
		}else{
			p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcol3D],y=UpRD$Aux[[1]][,input$ycol3D],z=UpRD$Aux[[1]][,input$zcol3D],color=UpRD$Aux[[1]][,input$catv],
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
		#withr::with_dir(file.path(UpRD$DivAna[[3]],"Output_3DPlots"),saveWidget(p,paste0('MDS3d_',input$catv3D,'.html'), selfcontained = F))
		p
	}
  })  
  
  output$try3d=renderPlotly({
	mds3Plot()
  })
##################################################################################################################################################################
#distance matrix heatmap
##################################################################################################################################################################
  heatPlot<-reactive({  
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	#req(mdata1())
	if(SelFile=="Data"){	  
		use=as.data.frame(BiodivInfo$res2[[1]])
		group=input$catv
		use=use[order(use[,group]),]
		useorder=match(use$Gen,rownames(BiodivInfo$res1[[6]]))
		distplot=plot_ly(x=rownames(BiodivInfo$res1[[6]])[useorder],y=rownames(BiodivInfo$res1[[6]])[useorder],z = BiodivInfo$res1[[6]][useorder,useorder], colorscale=input$colorheat,type = "heatmap")%>%
		layout(xaxis = list(showticklabels = F), yaxis = list(showticklabels = F))
		#withr::with_dir(file.path(BiodivInfo$res1[[3]]),saveWidget(distplot,'DistancesPlotEdited.html', selfcontained = F))
		distplot
	}else{
		UpRD=GenInfo$UploadRd
		use=as.data.frame(UpRD$Aux[[1]])
		group=input$catv
		use=use[order(use[,group]),]
		useorder=match(use$Gen,rownames(UpRD$DivAna[[6]]))
		distplot=plot_ly(x=rownames(UpRD$DivAna[[6]])[useorder],y=rownames(UpRD$DivAna[[6]])[useorder],z = UpRD$DivAna[[6]][useorder,useorder], colorscale=input$colorheat,type = "heatmap")%>%
		layout(xaxis = list(showticklabels = F), yaxis = list(showticklabels = F))
		#withr::with_dir(file.path(UpRD$DivAna[[3]]),saveWidget(distplot,'DistancesPlotEdited.html', selfcontained = F))
		distplot
	}
  })
  
  output$heat=renderPlotly({  
   heatPlot()
  })
  
##################################################################################################################################################################
#dendogram plot
################################################################################################################################################################## 
dendoPlot<-reactive({  
  if(!is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=1;SelFile="Data"}
	 if(is.null(GenInfo$dfgen) & !is.null(GenInfo$UploadRd)){Gdata=1;SelFile="RData"}
	 if(is.null(GenInfo$dfgen) & is.null(GenInfo$UploadRd)){Gdata=0}
	 hasResult <- !is.null(BiodivInfo$res1)
	 validate(
      need(Gdata != 0, "Please upload data"),
	  need(hasResult, "Please select options and click the button 'Run analysis'")
	)
	if(SelFile=="Data"){	
		#DivAna=BiodivInfo$res1
		#Aux=BiodivInfo$res2
		#save(DivAna,Aux,file="DivAna.RData")
		#if(!file.exists("Output_Dendograms")) dir.create("Output_Dendograms")
		data=as.data.frame(BiodivInfo$res2[[1]])
		info<- data[,c("Gen",input$catv)]
		info<- cbind(ID=info$Gen,info)
		names(info)=c("ID","Gen","Group")		
		tree=BiodivInfo$res2[[2]]
		
		if(!is.null(BiodivInfo$res1[[9]])){		
			q_df <- as.data.frame(BiodivInfo$res1[[9]])
			q_df$Individual <- rownames(q_df)
			q_long <- pivot_longer(q_df,
                cols = -Individual,
                names_to = "Cluster",
                values_to = "Mixture")						
		}
		
		if(typeof(info$Group)!="double"){
			if(!is.null(BiodivInfo$res1[[9]])){
                    updateSelectInput(session,'typeclust','Type',choices=c('circular'), selected='circular')			
					colp=input$colordend
				#save(tree,info,q_long,colp,file="infotree.RData")
					p <- ggtree(tree, layout = "circular")%<+%  info +
					geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, fontface="bold", show.legend=FALSE) +
					geom_point(aes(color = Group), alpha=0, show.legend=T) +
					scale_color_manual(values = input$colordend, breaks=levels(info$Group)) +
					geom_fruit(data = q_long,geom = geom_bar,mapping = aes(y = Individual,x = Mixture,fill = Cluster),
					orientation = "y",stat = "identity",width = 0.8,offset=0.01)+				
					theme(legend.position = input$poslen,
					legend.title = element_text(size = 12, face = "bold"),
					legend.text = element_text(size = 10),
					legend.key.size = unit(0.6, "cm")
					) + 
					guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1),ncol=2),
					fill = guide_legend(ncol = 2)) + 
					labs(color="Group", fill="Mixture") 
					p + xlim(0, max(p$data$x) + 0.5)				
			}else{
				p<-ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
					geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, fontface="bold", show.legend=FALSE) +
					geom_point(aes(color = Group), alpha=0, show.legend=T) +
					scale_color_manual(values=input$colordend, breaks=levels(info$Group))+					
					theme(legend.position = input$poslen,
					legend.title = element_text(size = 12, face = "bold"),
					legend.text = element_text(size = 10),
					legend.key.size = unit(0.6, "cm")
					) +
					guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1),ncol=2))+
					labs(color="Group") 
					p + xlim(0, max(p$data$x) + 0.5)	
			}
			
		}else{
			p<-ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
					scale_color_gradientn(colours=c("blue","cyan","green","orange","red")) +
					geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, fontface="bold", show.legend=FALSE) +
					geom_point(aes(color = Group), alpha=0, show.legend=T) +
					theme(legend.position = input$poslen,
					legend.title = element_text(size = 12, face = "bold"),
					legend.text = element_text(size = 10),
					legend.key.size = unit(0.6, "cm")
					) +
					guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1),ncol=2))+
					labs(color="Group") 
					p + xlim(0, max(p$data$x) + 0.5)	
		}		
		
	}else{
		UpRD=GenInfo$UploadRd
		#if(!file.exists("Output_Dendograms")) dir.create("Output_Dendograms")
		data=as.data.frame(UpRD$Aux[[1]][1])
		info<- data[,c("Gen",input$catv)]
		info<- cbind(ID=info$Gen,info)
		names(info)=c("ID","Gen","Group")
		kbest<-length(unique(info$Group))
		tree=UpRD$Aux[[2]]
		if(!is.null(UpRD$Aux[[1]][9])){
			q_df <- as.data.frame(UpRD$Aux[[1]][9])
			q_df$Individual <- rownames(q_df)
			q_long <- pivot_longer(q_df,
                cols = -Individual,
                names_to = "Cluster",
                values_to = "Mixture")			
		}
		
		if(typeof(info$Group)!="double"){
			if(!is.null(UpRD$Aux[[1]][9])){	
				   #updateSelectInput(session,'typeclust','Type',choices=c('circular','rectangular'), selected='circular')			
					p <- ggtree(tree, layout = "circular")%<+% info +
					geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, fontface="bold", show.legend=FALSE) +
					geom_point(aes(color = Group), alpha=0, show.legend=T) +
					scale_color_manual(values = input$colordend, breaks=levels(info$Group)) +
					geom_fruit(data = q_long,geom = geom_bar,mapping = aes(y = Individual,x = Mixture,fill = Cluster),
					orientation = "y",stat = "identity",width = 0.8,offset=0.01)+				
					theme(legend.position = input$poslen,
					legend.title = element_text(size = 12, face = "bold"),
					legend.text = element_text(size = 10),
					legend.key.size = unit(0.6, "cm")
					) + 
					guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1),ncol=2),
					fill = guide_legend(ncol = 2)) + 
					labs(color="Group", fill="Mixture") 
					p + xlim(0, max(p$data$x) + 0.5)				
			}else{
				p<-ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
					scale_color_manual(values=input$colordend)+
					geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, fontface="bold", show.legend=FALSE) +
					geom_point(aes(color = Group), alpha=0, show.legend=T) +
					theme(legend.position = input$poslen,
					legend.title = element_text(size = 12, face = "bold"),
					legend.text = element_text(size = 10),
					legend.key.size = unit(0.6, "cm")
					) +
					guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1),ncol=2))+
					labs(color="Group") 
					p + xlim(0, max(p$data$x) + 0.5)	
			}	
		}else{
			p<-ggtree(tree, layout=input$typeclust ,size=input$sizeline) %<+% info +
					scale_color_gradientn(colours=c("blue","cyan","green","orange","red")) +
					geom_tiplab(aes(label=Gen,color=Group),size=input$sizelab, offset=input$space, fontface="bold", show.legend=FALSE) +
					geom_point(aes(color = Group), alpha=0, show.legend=T) +
					theme(legend.position = input$poslen,
					legend.title = element_text(size = 12, face = "bold"),
					legend.text = element_text(size = 10),
					legend.key.size = unit(0.6, "cm")
					) +
					guides(color = guide_legend(override.aes = list(shape = 15, size = 5, alpha=1),ncol=2))+
					labs(color="Group") 
					p + xlim(0, max(p$data$x) + 0.5)	
		}
		
	}
	})
	output$dend=renderPlot({dendoPlot()})
  
    ## Report tab
    output$seeReportBio <- renderUI({
		req(report())
		tags$iframe(src = report(), width = "100%", height = "1100px", style = "border:none;")
	})
       
	 
	##Download distance matrix
    output$download_distM <- downloadHandler(
			filename = function() {
				paste0("DistanceMat_", Sys.Date(), ".csv")
			},
			content = function(file) {	
				req(BiodivInfo$res1)			
				mrdMAT1=BiodivInfo$res1[[6]]
				colnames(mrdMAT1)=seq(1:dim(BiodivInfo$res1[[6]])[1])			
				mrdMAT1=cbind(ID=seq(1:dim(BiodivInfo$res1[[6]])[1]),NAME=rownames(BiodivInfo$res1[[6]]),round(mrdMAT1,5))
				write.csv(mrdMAT1,file,quote=FALSE,row.names=FALSE)						
			}
	)
	
	##Download dendoplot
    output$download_dendPlot <- downloadHandler(
			filename = function() {
				paste0("DendogramPlot_", Sys.Date(), ".pdf")
			},
			content = function(file) {					
				p<-dendoPlot()
				req(p, p$data$x)
				p <- p + xlim(0, max(p$data$x) + 0.5)					
				ggsave(file,p,device="pdf")
			}
	)
	
	##Download report
    output$download_report <- downloadHandler(
		filename = function() {paste0("Biodiversity_dashboard_", Sys.Date(), ".html")},
  
		content = function(file) {
			req(report())
			file.copy(file.path("www", report()), file,overwrite = TRUE)
		}
	)
	 
	 report <- reactiveVal(NULL)
    
    observeEvent(input$generate_report, {

		shinybusy::show_modal_spinner(spin = "fading-circle", text = "Generating Report..." )

		src <- normalizePath("local/reportBio.Rmd")

		owd <- getwd()              # ⚠️ guardas ruta real
		setwd(tempdir())            # trabajas en temp

		file.copy(src, "reportBio.Rmd", overwrite = TRUE)

		outReport <- tryCatch({

			rmarkdown::render( "reportBio.Rmd",
				params = list(
					toDownload = TRUE,
					dend = dendoPlot(),
					mds2 = mds2Plot(),
					mds3 = mds3Plot(),
					heat = heatPlot(),
					amovaT = amovaTable(),
					popstr = DoforPopStr()[[1]],
					summBio = BiodivInfo
				),
				output_format = rmdformats::robobook(toc_depth = 4)
			)
		}, error = function(e) {
		showNotification(paste("Error:", e$message), type = "error")
		NULL
	})

  # 🔥 REGRESAR AL APP DIR
		setwd(owd)
		if (!is.null(outReport)) {
			if (!dir.exists("www")) dir.create("www")
			file.copy(outReport, "www/report.html", overwrite = TRUE)
			report("report.html")  # 🔥 clave
			shinyjs::click("download_report")
			shinyjs::click("download_distM")
			shinyjs::click("download_dendPlot")
		}
		shinybusy::remove_modal_spinner()
	})
    
	
  ########################################################################################################################################################################################    
  #Core-BIO-R
  ########################################################################################################################################################################################    
  #shinyFileChoose(input, 'filedistbio', roots = getVolumes(),filetypes=c('', 'csv'))
  #shinyFileChoose(input, 'filephendatbio', roots = getVolumes(),filetypes=c('', 'csv'))
  
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
  
  CoreInfo <- reactiveValues(genos1=NULL,statsF=NULL)

  datacore=reactive({
    req(input$datause)
    if ("phendat"%in%input$datause){
      dirfilePhen<-input$filephendatbio
    }else{      
      dirfilePhen<-"none"      
    }
    
    if ("gendat"%in%input$datause){
      datosgen=as.data.frame(GenInfo$dfgen)
      dirfileGen<-"exist"
      typedata<- cambia_caracter(input$typedata)	  
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
      dirfileDist<-input$filedistbio
    }else{
      dirfileDist<-"none" 
    }
    
    ver=sum(as.numeric(c(input$mrdEN,input$csedEN,input$gdEN,input$pcdEN,
                         input$mrdEE,input$csedEE,input$gdEE,input$pcdEE,
                         input$mrdAN,input$csedAN,input$gdAN,input$pcdAN,
                         input$SH,input$HE,input$CV)))
    
	### Correr funciones ----------------------------------------
	withProgress(message = 'Getting...', value = 0,{
	incProgress(1/2, detail = "Wait, Please!")
	
    resCH<-funchunt(datos,dirfileGen,dirfilePhen,dirfileDist,input$score,input$mrdEN,input$csedEN,input$gdEN,
             input$pcdEN,input$mrdAN,input$csedAN,input$gdAN,input$pcdAN,input$mrdEE,input$csedEE,input$gdEE,input$pcdEE,
             input$SH,input$HE,input$CV)
    	
	genos<-as.data.frame(resCH[[1]])
	distMat<-resCH[[2]]
	if (dirfileGen[1]!="none" & dirfilePhen[1]=="none" & dirfileDist[1]=="none"){
		rownames(datos)<-datos[,1]
		datos<-datos[,-1]
		statsF<-data.frame(cbind(c("He","sd_He","Ho","sd_Ho","CV","SH_Index","sd_SH"),rep(c("Pop_Complete","Core_Selected"),each=7),c(report_stats(datos),report_stats(datos[which(rownames(datos)%in%genos[,1]==T),]))))
		colnames(statsF)<-c("Parameter","Environment","Value")
	}

	if (dirfileGen[1]!="none" & dirfilePhen[1]!="none" & dirfileDist[1]=="none"){
		rownames(datos)<-datos[,1]
		datos<-datos[,-1]
		statsF<-data.frame(cbind(c("He","sd_He","Ho","sd_Ho","CV","SH_Index","sd_SH"),rep(c("Pop_Complete","Core_Selected"),each=7),c(report_stats(datos),report_stats(datos[which(rownames(datos)%in%genos[,1]==T),]))))
		colnames(statsF)<-c("Parameter","Environment","Value")
	}

	if (dirfileGen[1]!="none" & dirfilePhen[1]!="none" & dirfileDist[1]!="none"){
		rownames(datos)<-datos[,1]
		datos<-datos[,-1]
		statsF<-data.frame(cbind(c("He","sd_He","Ho","sd_Ho","CV","SH_Index","sd_SH"),rep(c("Pop_Complete","Core_Selected"),each=7),c(report_stats(datos),report_stats(datos[which(rownames(datos)%in%genos[,1]==T),]))))
		colnames(statsF)<-c("Parameter","Environment","Value")
	}

	if (dirfileGen[1]=="none" & dirfilePhen[1]!="none" & dirfileDist[1]=="none"){
		statsF<-data.frame(matrix(rep("Info no available",9),3,3))
		colnames(statsF)<-c("Parameter","Environment","Value")
	}

	if (dirfileGen[1]=="none" & dirfilePhen[1]=="none" & dirfileDist[1]!="none"){
		statsF<-data.frame(matrix(rep("Info no available",9),3,3))
		colnames(statsF)<-c("Parameter","Environment","Value")
	}
	
	todos <- unique(c(genos[,1], rownames(distMat)))
	t2 <- ifelse(todos %in% genos[,1] & todos %in% rownames(distMat), "selected", "All")
	genos1<-data.frame(gen = todos, Group = t2)

	##graphical representation, the MDS (multidimensional scaling) analysis
	mds=cmdscale(distMat, k=3,eig=T)
	coord=mds$points
	perctCP12=c(round(mds$eig[1]/sum(mds$eig)*100,2),round(mds$eig[2]/sum(mds$eig)*100,2),round(mds$eig[3]/sum(mds$eig)*100,2))
	colnames(coord)=c("Factor1","Factor2","Factor3")
	rm(mds)
		
	coord=cbind(gen=rownames(coord),coord)
	genos1<-merge(coord,genos1,"gen")
	genos1$gen=as.factor(genos1$gen)
	genos1$Factor1=as.numeric(genos1$Factor1)
	genos1$Factor2=as.numeric(genos1$Factor2)
	genos1$Factor3=as.numeric(genos1$Factor3)
	genos1$Group=as.factor(genos1$Group)
	incProgress(1, detail = "Finish")
	Sys.sleep(1)
	})
	
	vars=names(genos1)
    #Actualiza el selectinput de acuerdo a la base de datos cargada
    updateSelectInput(session,'xcolCH', 'X Variable',choices = vars,selected=vars[2])
    updateSelectInput(session,'ycolCH', 'Y Variable',choices = vars,selected=vars[3])
    updateSelectInput(session,'zcolCH', 'Z Variable',choices = vars,selected=vars[4])
    updateSelectInput(session,'catvCH', 'Group',choices = vars,selected=vars[5])
    updateSelectInput(session,'etiCH', 'Label',choices = vars,selected=vars[1])
	updateTextInput(session,'txCH','X Axis Label',value = paste0('Factor 1 (',perctCP12[1],'%)'))
	updateTextInput(session,'tyCH','Y Axis Label',value = paste0('Factor 2 (',perctCP12[2],'%)'))
	updateTextInput(session,'tzCH','Z Axis Label',value = paste0('Factor 3 (',perctCP12[3],'%)'))			
		
	#return(list(genos1,statsF))
	CoreInfo$genos1=genos1
	CoreInfo$statsF=statsF
    
  })
      
	output$defaultcore <- DT::renderDataTable({	
		validate(      
			need(!is.null(CoreInfo$statsF),"Select parameters and then click in Generate report button")
		)  	
		seedf<-as.data.frame(CoreInfo$statsF)
		datatable(seedf, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 7,width="100%", scrollX = TRUE))
	})

    output$selInd <- DT::renderDataTable({   
		validate(      
			need(!is.null(CoreInfo$genos1),"Select parameters and then click in Generate report button")
		)
		seedf<-as.data.frame(CoreInfo$genos1)
		datatable(seedf, selection="multiple", escape=FALSE, 
              options = list(sDom  = '<"top">lrt<"bottom">ip',pageLength = 10,width="100%", scrollX = TRUE))
	})
	
observeEvent(input$catvCH,{	
    genos1<-CoreInfo$genos1
	set.seed(7)
    colores=colors()[-c(1,3:12,13:25,24,37:46,57:67,80,82,83,85:89,101:106,108:113,126:127,138,140:141,152:253,260:366,377:392,
                        394:447,449,478:489,492,513:534,536:546,557:561,579:583,589:609,620:629,418,436,646:651)]
    
	if(typeof(genos1[,input$catvCH])!="double"){
		d=sample(colores,100)
		var=as.factor(genos1[,input$catvCH])
		grupos=nlevels(var)
	}else{
		d=c("Jet","Jet","Jet","Jet","Jet","Jet")
		grupos=2
	}
	
    updateSelectInput(session,'colorCH','Choose a color', choices=d,selected=d[1:grupos])
  })
##################################################################################################################################################################
  #grafico 2d Core Hunter
##################################################################################################################################################################
  mds2PlotCore<-reactive({
		validate(      
			need(!is.null(CoreInfo$genos1),"Select parameters and then click in Generate report button")
		)
		genos1<-CoreInfo$genos1
		if(typeof(genos1[,input$catvCH])!="double"){
			p=plot_ly(data=genos1,x=genos1[,input$xcolCH],y=genos1[,input$ycolCH],color=genos1[,input$catvCH],
					type="scatter",mode="markers",colors = input$colorCH,xaxis=F, yaxis=F,
					text=genos1[,input$etiCH],marker=list(size=input$sizeCH))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgpCH)%>%
			#titulo y etiquetas ejes
			layout(title=input$tpCH,titlefont=list(size=input$tsCH,color=input$pncCH), xaxis = list(title = input$txCH, titlefont=list(size=input$szlCH,color=input$acCH)),
					yaxis = list(title = input$tyCH,titlefont=list(size=input$szlCH,color=input$acCH)))
		}else{
			p=plot_ly(data=genos1,x=genos1[,input$xcolCH],y=genos1[,input$ycolCH],color=genos1[,input$catvCH],
					type="scatter",mode="markers",xaxis=F, yaxis=F,colors=c("blue","cyan","green","orange","red"),
					text=genos1[,input$etiCH],marker=list(size=input$sizeCH))%>%
			#color de fondo del grafico
			layout(plot_bgcolor=input$bkgpCH)%>%
			#titulo y etiquetas ejes
			layout(title=input$tpCH,titlefont=list(size=input$tsCH,color=input$pncCH), xaxis = list(title = input$txCH, titlefont=list(size=input$szlCH,color=input$acCH)),
					yaxis = list(title = input$tyCH,titlefont=list(size=input$szlCH,color=input$acCH)))		
		}		
		p
	#}else{
	#	UpRD=GenInfo$UploadRd
	#	if(!file.exists("Output_2DPlots")) dir.create("Output_2DPlots")    
	#	if(typeof(UpRD$Aux[[1]][,input$catvCH])!="double"){
	#		p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcolCH],y=UpRD$Aux[[1]][,input$ycolCH],color=UpRD$Aux[[1]][,input$catvCH],
	#				type="scatter",mode="markers",colors = input$colorCH,xaxis=F, yaxis=F,
	#				text=UpRD$Aux[[1]][,input$etiCH],marker=list(size=input$sizeCH))%>%
	#		#color de fondo del grafico
	#		layout(plot_bgcolor=input$bkgpCH)%>%
	#		#titulo y etiquetas ejes
	#		layout(title=input$tpCH,titlefont=list(size=input$tsCH,color=input$pncCH), xaxis = list(title = input$txCH, titlefont=list(size=input$szlCH,color=input$acCH)),
	#				yaxis = list(title = input$tyCH,titlefont=list(size=input$szlCH,color=input$acCH)))
	#	}else{
	#		p=plot_ly(data=UpRD$Aux[[1]],x=UpRD$Aux[[1]][,input$xcolCH],y=UpRD$Aux[[1]][,input$ycolCH],color=UpRD$Aux[[1]][,input$catvCH],
	#				type="scatter",mode="markers",xaxis=F, yaxis=F,colors=c("blue","cyan","green","orange","red"),
	#				text=UpRD$Aux[[1]][,input$etiCH],marker=list(size=input$sizeCH))%>%
	#		#color de fondo del grafico
	#		layout(plot_bgcolor=input$bkgpCH)%>%
	#		#titulo y etiquetas ejes
	#		layout(title=input$tpCH,titlefont=list(size=input$tsCH,color=input$pncCH), xaxis = list(title = input$txCH, titlefont=list(size=input$szlCH,color=input$acCH)),
	#				yaxis = list(title = input$tyCH,titlefont=list(size=input$szlCH,color=input$acCH)))
	#	}
	#	#el siguiente codigo, cambia temporalmente el directorio de trabajo para guardar el grafico 2d
	#	#primero se especifica la direccion en la que se guardara y luego la accion (guardar el grafico)
	#	withr::with_dir(file.path(UpRD$DivAna[[3]],"Output_2DPlots"),saveWidget(p,paste0('MDS2dCH_',input$catvCH,'.html'), selfcontained = F))
	#	p
	#}
  })
  
  output$tryCH=renderPlotly({   
	mds2PlotCore()
  })
  
  valdatacore <- reactive({    
    req(input$datause)
    # Ejemplo simple
    dirfilePhen <- input$filephendatbio
    dirfileDist <- input$filedistbio
    
    weights <- as.numeric(c(
      input$mrdEN,input$csedEN,input$gdEN,input$pcdEN,
      input$mrdEE,input$csedEE,input$gdEE,input$pcdEE,
      input$mrdAN,input$csedAN,input$gdAN,input$pcdAN,
      input$SH,input$HE,input$CV
    ))
    
    ver <- sum(weights, na.rm = TRUE)
    # NO validate aquí
    return(list(
      ver = ver,
      dirfilePhen = dirfilePhen,
      dirfileDist = dirfileDist
    ))
  })
  
  
observeEvent(input$generate_core, {

  req(valdatacore())
  data <- valdatacore()

  # 🔍 Validaciones (tu lógica)
  errores <- c()

  if(("phendat" %in% input$datause) && is.null(data$dirfilePhen)){
    errores <- c(errores, "Please select phenotypic data")
  }

  if(("distdat" %in% input$datause) && is.null(data$dirfileDist)){
    errores <- c(errores, "Please select distance data")
  }

  if(("gendat" %in% input$datause) && is.null(GenInfo$dfgen)){
    errores <- c(errores, "Please select genetic data")
  }

  if(!("gendat" %in% input$datause) && ("distdat" %in% input$datause)){
    errores <- c(errores, "Combination of data NOT available")
  }

  if(data$ver != 1){
    errores <- c(errores, paste("The sum of weights must be 1 (actual:", data$ver, ")"))
  }

  # 🚫 Si hay errores → no generar
  if(length(errores) > 0){
    #showNotification(paste(errores, collapse = "\n"), type = "warning")
	showModal(modalDialog(title = "Warnings", paste(errores, collapse = "\n"),    
    easyClose = TRUE,
    footer = modalButton("Close"),
    size = "l"  # 👈 tamaño grande
    ))
    return()
  }

  # 🔥 Generar reporte
  shinybusy::show_modal_spinner(text = "Generating Report...")

		datacore()

	shinybusy::remove_modal_spinner()
	
  })
  
  ##Download report
    output$download_core <- downloadHandler(
		filename = function() {paste0("CoreSubset_dashboard_", Sys.Date(), ".html")},
  
		content = function(file) {
			req(reportCore())
			file.copy(file.path("www", reportCore()), file,overwrite = TRUE)
		}
	)
	 
	 reportCore <- reactiveVal(NULL)
    
    observeEvent(input$downCore,{
		
		shinybusy::show_modal_spinner(spin = "fading-circle", text = "Download Report..." )
		
		src <- normalizePath("local/reportCore.Rmd")

		owd <- getwd()              # ⚠️ guardas ruta real
		setwd(tempdir())            # trabajas en temp

		file.copy(src, "reportCore.Rmd", overwrite = TRUE)

		outReport <- tryCatch({

			rmarkdown::render( "reportCore.Rmd",
				params = list(
					toDownload = TRUE,
					resCore = list(genos1=CoreInfo$genos1,statsF=CoreInfo$statsF),
					mds2 = mds2PlotCore()					
				),
				output_format = rmdformats::robobook(toc_depth = 4)
			)
		}, error = function(e) {
		showNotification(paste("Error:", e$message), type = "error")
		NULL
	})

  # 🔥 REGRESAR AL APP DIR
		setwd(owd)
		if (!is.null(outReport)) {
			if (!dir.exists("www")) dir.create("www")
			file.copy(outReport, "www/reportCore.html", overwrite = TRUE)
			reportCore("reportCore.html")  # 🔥 clave	
			shinyjs::click("download_core")			
		}
		shinybusy::remove_modal_spinner()
	
})
    

  
  
shinyjs::disable("downCore")

observe({
  if (!is.null(CoreInfo$genos1)) {
    shinyjs::enable("downCore")
  }
})
  
##################################################################################################################################################################
###################################################################################################################################################################################################
#Start GWAS Code
##################################################################################################################################################################
##################################################################################################################################################################
#shinyFileChoose(input, 'filephen', roots = getVolumes(),filetypes=c('csv'))
#shinyFileChoose(input, 'fileplants', roots = getVolumes(),filetypes=c('gz'))

myDataPhen<-reactive({
	req(GenInfo)
	inFilephen= input$filephen
	if(is.null(inFilephen)){PhenData=0}else{PhenData=1}	 
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
	req(GenInfo)
	inFileplants=input$fileplants
	if(is.null(inFileplants)){PlantsData=0}else{PlantsData=1}
	validate(
		need(PlantsData != 0, "Please select genome ref data")
    )
	print("Please wait...")
	#dirdata=str_replace(input$filegen$datapath,input$filegen$name,"")
	#load(paste0(dirdata,"Output_BIO-R\\ChromMap_",cambia_caracter(strsplit(input$filegen$name,"[.]","")[[1]][1]),"\\MapChrom.RData"))	
	temp_file <- file.path(tempdir(), "archivotmp.gtf.gz")
	file.copy(input$fileplants$datapath, temp_file, overwrite = TRUE)

	wheat<-import(temp_file)
	wheat=as.data.frame(wheat)	
	check=sort(as.character(unique(wheat[,1])[1:dim(GenInfo$chromosome_file)[1]]))
	checknum=as.numeric(levels(as.factor(check)))
	wheat[,1]=as.character(wheat[,1])
	for (k in 1:dim(GenInfo$chromosome_file)[1]){
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
		#dirdata=str_replace(input$filegen$datapath,input$filegen$name,"")
		#load(paste0(dirdata,"Output_BIO-R\\ChromMap_",cambia_caracter(strsplit(input$filegen$name,"[.]","")[[1]][1]),"\\MapChrom.RData"))
		
		datos<-GenInfo$hapmap
		datos<-datos[,-c(2:11)]
		newcolnames <- cambia_caracter(quita_espacio(as.matrix(colnames(datos))))
		colnames(datos) <- putg(newcolnames)
		datos<-datos[which(datos[,1]%in%GenInfo$annotation_file$ID==T),]
		rownames(datos)<-datos[,1]
		datos<-t(datos[,-1])
		
		posit<-GenInfo$posit
		posit<-posit[,c(3,1,2)]
		posit<-posit[which(posit[,1]%in%GenInfo$annotation_file$ID==T),]
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
		
		#outFolderPh <- cambia_caracter(paste0("GWASAnalysis_",strsplit(input$filegen$name,"[.]","")[[1]][1]))
		#setwd(dirdata)
		#if(!file.exists("Output_BIO-R")) dir.create("Output_BIO-R")
		#setwd("Output_BIO-R")
		#if(!file.exists(outFolderPh)) dir.create(outFolderPh)
		#setwd(outFolderPh)
		
		#save(gDataDrops,GWASDrops,file="DataGWAS.RData")
		
		#for(nt in 1:length(Ttraits)){			
		#	ass=split(GWASDrops[["GWAResult"]][["dropsPheno"]],GWASDrops[["GWAResult"]][["dropsPheno"]][["trait"]])
		#	assL[[nt]]=ass[Ttraits[nt]]
			#write.csv(ass,paste0("MarkersEffect_",Ttraits[nt],".csv"))		
		#}
	return(list(GWASDrops))
	}
})
	

DoforGene<-reactive({
	if(input$startAna=="StarChrom" & input$typedata=="vcfile"){
	req(myDataRegGen())
	wheat<-as.data.frame(myDataRegGen()[[1]])
	#chrdim<-as.data.frame(myDataRegGen()[[2]])
	GWASDrops<-DoforGWAS()[[1]]	
	#save(GWASDrops,wheat,file="DataGWASRegGen.RData")
	
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
	
	#savein<-DoforGWAS()[[2]]
	#setwd(savein)	
	
	return(list(ass,wheat,htm))
	}
})

##################################################################################################################################################################
#For see regional associated gene plot
##################################################################################################################################################################  
  RegGenplotR<-reactive({
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
	}else{
		load(input$fileRData$datapath)
	}
	
	wheatmp=wheat[which(wheat[,1]==as.character(input$RegGenPlotChr)),]
	
	IntReg=IntRegionalPlot(chr=as.numeric(input$RegGenPlotChr),left=as.numeric(input$RegGenPlotleft),
	right=as.numeric(input$RegGenPlotrigth),threshold=as.numeric(input$RegGenPlotChrThr),gtf=wheat,association=ass,
	hapmap=htm,label_gene_name=TRUE,hapmap_ld=htm,leadsnp_size=3)
	listgen=IntReg[["plot_env"]][["gene_list"]]	
	
	return(IntReg,listgen)	
  })
  
  output$RegGenPlotChrPlot=renderPlot({
	IntReg=RegGenplotR()[[1]]
	print(IntReg)
  })
##################################################################################################################################################################
#For see manhattan plot
##################################################################################################################################################################  
  ManplotR<-reactive({
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
	}else{
		load(input$fileRData$datapath)
	}
	
	if(input$Manchr=="All"){ManchrV=1:length(unique(GWASDrops[["GWAResult"]][["dropsPheno"]][["chr"]]))}else{ManchrV=eval(parse(text=input$Manchr))}
	if(input$MancolPalette=="Blue"){MancolorPaletteV=rep(c("blue","darkblue"),length(ManchrV))[1:length(ManchrV)]}
	if(input$MancolPalette=="Gray"){MancolorPaletteV=rep(c("black","cornsilk4"),length(ManchrV))[1:length(ManchrV)]}
	if(input$MancolPalette=="Green"){MancolorPaletteV=rep(c("darkolivegreen1","darkolivegreen4"),length(ManchrV))[1:length(ManchrV)]}
	if(input$MancolPalette=="Colors"){MancolorPaletteV=colorRampPalette(brewer.pal(11,"Spectral"))(length(ManchrV))}
	
	p<-plot(GWASDrops, plotType = "manhattan", trait = input$ManTrait, yThr = as.numeric(input$ManyThr) ,chr= ManchrV ,colPalette= MancolorPaletteV)	
	p<-p + ggtitle(paste0("Trait: ", input$ManTrait,"  LOD: ", input$ManyThr))
	p				   	
  })
  
  output$ManPlot=renderPlot({
	ManplotR()
  })
##################################################################################################################################################################
#For see QTL plot
################################################################################################################################################################## 
  QTLplotR<-reactive({
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
	}else{
		load(input$fileRData$datapath)
	}
	
	if(input$QTLchr=="All"){QTLchrV=1:length(unique(GWASDrops[["GWAResult"]][["dropsPheno"]][["chr"]]))}else{QTLchrV=eval(parse(text=input$QTLchr))}
	
	p<-plot(GWASDrops, plotType = "qtl", yThr = as.numeric(input$QTLyThr) ,chr= QTLchrV)	
	p<-p + ggtitle(paste0("LOD: ", input$QTLyThr))
	p				   
  })
  
  output$QTLPlot=renderPlot({
	QTLplotR()
  })
##################################################################################################################################################################
#For see QQ plot
##################################################################################################################################################################  
  QQplotR<-reactive({  
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
	}else{
		load(input$fileRData$datapath)
	}
		
	#If the lambda value is greater than 1, then this may be evidence for some systematic bias that needs to be corrected in your analysis.
	print("Lambda for GWAS:")
	print(GWASDrops[["GWASInfo"]][["inflationFactor"]])
		
	p<-plot(GWASDrops, plotType = "qq", trait=input$QQTrait)
	p<-p + ggtitle(paste0("Trait: ", input$QQTrait, "   Lambda for GWAS: ", round(GWASDrops[["GWASInfo"]][["inflationFactor"]][[1]][input$QQTrait],4) ))
	p				   
})
	
	output$QQPlot=renderPlot({
		QQplotR()
	})


##Download report
    output$generateGWAS <- downloadHandler(
		filename = function() {paste0("GWAS_dashboard_", Sys.Date(), ".html")},
  
		content = function(file) {
			req(reportGWAS())
			file.copy(file.path("www", reportGWAS()), file,overwrite = TRUE)
		}
	)
	 
	 reportGWAS <- reactiveVal(NULL)
    
    observeEvent(input$generate_GWAS,{

		shinybusy::show_modal_spinner(spin = "fading-circle", text = "Generating Report..." )

		src <- normalizePath("local/reportGWAS.Rmd")

		owd <- getwd()              # ⚠️ guardas ruta real
		setwd(tempdir())            # trabajas en temp

		file.copy(src, "reportGWAS.Rmd", overwrite = TRUE)
		
		if(!is.null(input$fileplants$datapath)){

			outReport <- tryCatch({

				rmarkdown::render( "reportGWAS.Rmd",
					params = list(
						toDownload = TRUE,
						resGWAS = DoforGWAS(),
						qqplotR = QQplotR(),
						qtlplotR = QTLplotR(),
						manplotR = ManplotR(),
						rengenplotR = RenGenplotR()
					),
					output_format = rmdformats::robobook(toc_depth = 4)
				)
			}, error = function(e) {
				showNotification(paste("Error:", e$message), type = "error")
				NULL
			})
		} else{
			outReport <- tryCatch({

				rmarkdown::render( "reportGWAS.Rmd",
					params = list(
						toDownload = TRUE,
						resGWAS = DoforGWAS(),
						qqplotR = QQplotR(),
						qtlplotR = QTLplotR(),
						manplotR = ManplotR(),
						rengenplotR = NULL
					),
					output_format = rmdformats::robobook(toc_depth = 4)
				)
			}, error = function(e) {
				showNotification(paste("Error:", e$message), type = "error")
				NULL
			})
		}
	# 🔥 REGRESAR AL APP DIR
		setwd(owd)
		if (!is.null(outReport)) {
			if (!dir.exists("www")) dir.create("www")
			file.copy(outReport, "www/reportGWAS.html", overwrite = TRUE)
			reportGWAS("reportGWAS.html")  # 🔥 clave	
			shinyjs::click("generateGWAS")			
		}
		shinybusy::remove_modal_spinner()
	})

##################################################################################################################################################################
}



