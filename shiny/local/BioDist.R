#mrdMAT=read.csv("SNP_distance_matrix.csv")[,-1]
#colorMDS=read.csv("Groups10TID.csv")
#nclust=3
#NEWONE="ResultsPlots"

biodist<-function(mrdMAT,coloMDS,nclust){
NEWONE="DistMatRes"
if(!file.exists(NEWONE)) dir.create(NEWONE)
setwd(NEWONE)
ForPlotly=colorMDS
mrdMAT=as.matrix(mrdMAT)
rownames(mrdMAT)=putg(colnames(mrdMAT))
colnames(mrdMAT)=putg(colnames(mrdMAT))

suppressWarnings(library(Hmisc))  
suppressWarnings(library(plotly))  
suppressWarnings(library(cluster))
suppressWarnings(library(car))
suppressWarnings(library(ape))
suppressWarnings(library(RColorBrewer))


distplot=plot_ly(x=rownames(mrdMAT),y=rownames(mrdMAT),z = mrdMAT, colorscale="Reds",type = "heatmap")%>%
   layout(xaxis = list(showticklabels = F), yaxis = list(showticklabels = F),
       # color option buttons  
       updatemenus =list( list(
         type = "buttons",         
         buttons = list(           
           list(method = "restyle",
                args = list("colorscale", "Reds"),
                label = "Reds"),           
           list(method = "restyle",
                args = list("colorscale", "Jet"),
                label = "Jet"),           
           list(method = "restyle",
                args = list("colorscale", "Blues"),
                label = "Blues"),           
           list(method = "restyle",
                args = list("colorscale", "Viridis"),
                label = "Viridis")
         ))))
namedisplot=paste("DistancesPlot.html",sep="")
htmlwidgets::saveWidget(as.widget(distplot),selfcontained=FALSE,namedisplot)
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
clust=agnes(mrdMAT, method = "ward")
#########################################################################
#save(traits,mayorque,menorque,coord,clust,datos,exadivg,id1,mrdMAT,perctCP12,file=paste(dir_root,"/Programs/Rcode/tmpoutput.rda",sep=""))	

   ##save MDS file
   coord1=as.data.frame(coord)
   names(coord1)=c("Factor1","Factor2","Factor3")
   write.csv(coord1,"MDStable.csv",quote=FALSE)      
   ##Create file with the options that were chosen for the analysis

   nall=2
  out="SummaryDiversityAnalysis.csv"
   
   save.color=unique(rgb(t(col2rgb(colors())),max=255))
   save.color=save.color[c(2,3:15,23:28,38:57,62,63:97,103:112,114:122,124,127:135,235:253,263:306,308:348,355:392,394:436,442:451,461:469,474:493,499:502)]
   
  #if (colclust=="BOTH"){ 
       groupCodes=as.character(colorMDS[,2])
       colorCodes=save.color[sample(1:333,size=length(unique(groupCodes)),replace=FALSE)]		   
       names(colorCodes)=unique(groupCodes)	
       
       #########################################################################
       ##codigo para Archaeopteryx##############################################
       ######################################################################### 
       TFArx=as.phylo(as.hclust(clust))
       if (file.exists("proofnewick.txt")) file.remove("proofnewick.txt")  
       write.tree(TFArx,"proofnewick.txt")  
       line2=as.vector(scan("proofnewick.txt",what=character()))
       line1=substring(line2,1:nchar(line2),1:nchar(line2))
       groupCodesB=cutree(as.hclust(clust),k=nclust)
       colorCodesB=col2rgb(save.color[sample(1:333,size=nclust,replace=FALSE)])
       filearch=paste("TreeforArchaeopteryx.xml",sep="")           
       export.xml(filearch,line1,TFArx,2,colorCodes,groupCodes,colorCodesB,groupCodesB)
       filecopy="TreeforArchaeopteryx.xml"
       file.copy(filearch,filecopy)		   
       if (file.exists("proofnewick.txt")) file.remove("proofnewick.txt")
       #########################################################################	
     pp= as.data.frame(cutree (as.hclust(clust), k = nclust))
     groups=as.data.frame(pp)
     
     coord1=cbind(gen=rep(0,dim(coord)[1]),coord)
     
      coord1[,1]=groups[,1]
       for (i in 2:dim(ForPlotly)[2]){
         coord1=cbind(coord1,as.character(ForPlotly[,i]))
       }			
       namemds2D=paste("mds2D","BOTH.html",sep="")
       namemds3D=paste("mds3D","BOTH.html",sep="")
       coord1=as.data.frame(coord1)
       names(coord1)=c("Group","Factor1","Factor2","Factor3",names(ForPlotly[-1]))
       coord1$Group=as.factor(coord1$Group)
       coord1$Factor1=as.numeric(as.character(coord1$Factor1))
       coord1$Factor2=as.numeric(as.character(coord1$Factor2))
       coord1$Factor3=as.numeric(as.character(coord1$Factor3))
       numlevel=list()
       namesgroups=list()
       numlevel[[1]]=nlevels(coord1[,1])
       namesgroups[[1]]="Cluster"
       for (i in 1:(dim(ForPlotly)[2]-1)){
         coord1[,4+i]=as.factor(coord1[,4+i])
         numlevel[[i+1]]=nlevels(coord1[,4+i])
         namesgroups[[i+1]]=names(coord1)[4+i]
       }
       ## graphic representation for MDS two dimension
       f <- list(family = "Brodway", size = 18)
       x <- list(title = paste("CP1(",perctCP12[1],"%)",sep=""),titlefont = f)
       y <- list(title = paste("CP2(",perctCP12[2],"%)",sep=""),titlefont = f)
       z <- list(title = paste("CP3(",perctCP12[3],"%)",sep=""),titlefont = f)
       pString<-paste("mds2<-plot_ly(coord1,x = ~Factor1)%>%add_annotations(x = 0.47, y = 1.07, text = 'Points Color',xref = 'paper', yref = 'paper',showarrow = FALSE,font = list(color = 'red',size = 14))%>%  add_annotations(x = 0.37, y = 1.07, text = 'Groups',xref = 'paper', yref = 'paper',showarrow = FALSE,font = list(color = 'red',size = 14))%>%add_markers(y = ~Factor2,color=~",eval(names(coord1)[1]),",colors='Set1',mode='markers',marker=list(size=13),text= ~paste('Genotype: ', as.character(rownames(coord1))),visible=T)",sep="")
       
       for (i in 2:length(namesgroups)){			    
         pString<-paste(pString,"%>%add_markers(y = ~Factor2,color=~",eval(namesgroups[[i]]),",colors='Set1',mode='markers',marker=list(size=13),text= ~paste('Genotype: ', as.character(rownames(coord1))),visible=F)",sep="")
       }
       eval(parse(text=pString))
       allbuttons=list()
       vectFT=rep(TRUE,numlevel[1])
       for (i in 2:length(numlevel)){vectFT=c(vectFT,rep(FALSE,numlevel[[i]]))}			
       allbuttons[[1]]=list(method="restyle", 
                            args=list("visible",as.list(vectFT)),
                            label=namesgroups[[1]])			
       for (i in 2:length(numlevel)){		
         vectFT=FALSE			
         for (j in 1:length(numlevel)){
           if(i==j){
             vectFT=c(vectFT,rep(TRUE,numlevel[[j]]))
           }else{vectFT=c(vectFT,rep(FALSE,numlevel[[j]]))}
         }
         allbuttons[[i]]=list(method="restyle", 
                              args=list("visible",as.list(vectFT[-1])),
                              label=namesgroups[[i]])
       }			
       mds2<-mds2%>%layout(xaxis = x, yaxis = y, 
                           updatemenus=list( 
                             list(
                               active=-1,
                               type="buttons",  
                               direction = "right",
                               yanchor = "top",							
                               x = 0.4,
                               y = 1.052,							
                               buttons=allbuttons),
                             list(
                               active=-1,
                               type="buttons",
                               direction = "right",
                               yanchor = "top",														
                               x = 0.5,
                               y = 1.05,							
                               buttons=list(
                                 #Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn  Accent 8 Dark2 8 Paired 12 Pastel1 9 Pastel2 8 Set1 9 Set2 8 Set3 12 
                                 list(method = "restyle",
                                      args = list("marker.color", colorRampPalette(brewer.pal(10,"Spectral"))(max(unlist(numlevel)))),
                                      label = "SetColor1")
                               ))
                           ))                  
       
       htmlwidgets::saveWidget(as.widget(mds2),selfcontained=FALSE,namemds2D)
       
     
       coordcurly=cbind(paste("category",groups[,1],sep=" "),rownames(coord),coord1[,2:4])
       colnames(coordcurly)=c("categories:Cluster","Label","Factor1","Factor2","Factor3")	   
       write.table(coordcurly,"curlyfile.txt",sep = "\t" ,row.names=FALSE,quote=FALSE)
#######################################################################################################
#######################################################################################################
       
}