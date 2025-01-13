#Procedure amova
forAMOVA<-function(agc.dist,agc.env){
     #suppressWarnings(library(vegan))	 
     agc.adonis <- adonis2(agc.dist ~ Pop, data=agc.env, permutations=999)
     adonis.data <- data.frame(agc.adonis)
     adonis.data$MS <- adonis.data$SumOfSqs / adonis.data$Df     
     pop.sizes <- table(agc.env$Pop)
     n1 <- mean(pop.sizes)     
     # how much variance is detected at each stratification.
     #Variations  Within samples
     Sigma.withinSamples <- adonis.data["Residual", "MS"]     
     #Variations  Between samples
     Sigma.betweenSamples <- (adonis.data[1, "MS"] - Sigma.withinSamples) / n1     
     Sigma.total <- Sigma.withinSamples + Sigma.betweenSamples      
     #Phi-pop-total
     PHI.ST <- Sigma.betweenSamples / Sigma.total     
     variance=cbind(adonis.data[,c(1,2,6,4,5)],
                    Sigma=c(Sigma.betweenSamples,Sigma.withinSamples,Sigma.total),
                    PercVar=c(100*(Sigma.betweenSamples/Sigma.total),100*(Sigma.withinSamples/Sigma.total),100),
                    Phi=c(PHI.ST,NA,NA)
                    )
     variance=cbind(c("Between Groups","Within Groups","Total"),variance)
	 colnames(variance)[1]=c("source")
	 return(variance)
}

gdiv<-function(datos,datos1,selv,quitomono,mdrMAT){

######################################################
##Diversity among groups
######################################################
out="PopulationStructure.csv"
## function for the HE calculus
heiter=function(freqs) {  
    if(is.matrix(freqs)){
        f1=t(cbind(freqs[,1:ncol(freqs)],1-freqs[,1:ncol(freqs)]))
        nf1=t(t(apply(f1,1,function(x) sum(x,na.rm=T))))
        nf11=nf1[1:dim(freqs)[2]]
        nf12=nf1[dim(freqs)[2]+1:dim(freqs)[2]*2]
        suma1=nf11+nf12
        h0=1-((nf11/suma1)^2+(nf12/suma1)^2)
        he=mean(h0,na.rm=T) 
    }else{he=NA}
} 

#For calculate diversity
pp=as.data.frame(as.character(datos1[,selv]))	
rownames(pp)=datos1[,1]       
groups=pp        
ngroup=nlevels(as.factor(groups[,1]))
groups[,1]=as.numeric(as.factor(groups[,1]))

######################################################
##Diversity among groups by cluster
###################################################### 
       if(selv=="GroupClust"){
	        ## end of HE function
            HE=mean(1-((datos$pest)^2+(datos$pest1)^2))			
            tfreq=t(datos[,1:(ncol(datos)-4)])																	
            rownames(tfreq)=colnames(datos[,1:(ncol(datos)-4)])            
            nummark=list()
            heW=list()
            for (i in 1:ngroup){
              tfreq2=datos[,which(groups[,1]==i)]
              nacc=dim(tfreq2)[2]
              nmark=dim(tfreq2)[1]
			  if(length(which(groups[,1]==i))==1 || length(which(groups[,1]==i))==0){
                 nummark[[i]]=0
                 heW[[i]]="There is no group"
              }else{
                    tfreq2$pest=apply(tfreq2[,1:nacc],1,mean,na.rm=T)     
                    tfreq2$pest1=apply((1-tfreq2[,1:nacc]),1,mean,na.rm=T)
                    if (quitomono==TRUE){
					     #solo quito monomorficos
                         vect2=which(tfreq2$pest!=1 & tfreq2$pest1!=1 & tfreq2$pest!=0 & tfreq2$pest1!=0)
                         if (length(vect2)!=0){
                           tfreq2=tfreq2[vect2,]
                           nummark[[i]]=length(vect2)
                           tfreq2=t(tfreq2[,1:(ncol(tfreq2)-2)])
                           heW[[i]]=heiter(tfreq2)
                         }else{
                              nummark[[i]]=0
                              heW[[i]]="All markers monomorphics"
                              }
					}else{
					    vect2=which(tfreq2$pest!=1 & tfreq2$pest1!=1 & tfreq2$pest!=0 & tfreq2$pest1!=0)
						print(length(vect2))
						if (length(vect2)!=0){
						    writemono=cbind(rownames(tfreq2)[-vect2],tfreq2$pest[-vect2])
						    write.csv(writemono,paste0(getwd(),"/Output_MarkMonoGroups/MarkersMonomorphicsForGroup_",selv,"_in level_",levels(pp[,1])[i],".csv"),row.names=F)
						}
                        tfreq2=t(tfreq2[,1:(ncol(tfreq2)-2)])
						nummark[[i]]=nmark
                        heW[[i]]=heiter(tfreq2)  
					}
				}
            }                   
            heT=HE 																	## Population He              			
            heW=unlist(heW)
            nummark=c(dim(tfreq)[2],unlist(nummark))
            heB=heT-mean(heW,na.rm=TRUE) 																## Among individuals whitin groups He 
            Fst=heB/heT 																	## Diversity among groups (estadistico F wright el grado de 
                                                         ## diferencia genetica entre las poblaciones, en funcion de las frecuencias alelicas            
            div2=cbind(t(t(c("Diversity among groups","Diversity within group",rep("",(length(heW)-1))))),t(t(c("",seq(1:ngroup)))),t(t(c(Fst,heW))))
            div2=cbind(div2,nummark)
			
			agc.env=as.data.frame(groups[,1])
			names(agc.env)<-c("Pop")
			agc.env$Pop<-as.factor(agc.env$Pop)
			tabamv=forAMOVA(mdrMAT,agc.env)
            if("div2"%in%ls()==TRUE){
            cat("\n","Population Structure for group: ",selv,"\n","\n",file=out,append=T)
            write.table(div2, file = out, append = T,quote=F, sep=",",col.names=F,row.names=F)
			cat("\n","AMOVA for group: ",selv,"\n","\n",file=out,append=T)
            write.table(tabamv, file = out, append = T,quote=F, sep=",",col.names=T,row.names=F)			
            }
	   }
	   if(selv!="GroupClust"){                       			
	        ## end of HE function
            HE=mean(1-((datos$pest)^2+(datos$pest1)^2))
            tfreq=t(datos[,1:(ncol(datos)-4)])																	
            rownames(tfreq)=colnames(datos[,1:(ncol(datos)-4)])            
            nummark=list()
            heW=list()
			outheT=list()
            for (i in 1:ngroup){
              tfreq2=datos[,which(groups[,1]==i)]
              nacc=dim(tfreq2)[2]
              nmark=dim(tfreq2)[1]
			  if(length(which(groups[,1]==i))==0){
                 nummark[[i]]=0
                 heW[[i]]="There is no group"
			  }
			  if(length(which(groups[,1]==i))==1){
                 nummark[[i]]=0
                 heW[[i]]="There is no group"
				 outheT[[i]]=which(groups[,1]==i)
              }else{
                    tfreq2$pest=apply(tfreq2[,1:nacc],1,mean,na.rm=T)     
                    tfreq2$pest1=apply((1-tfreq2[,1:nacc]),1,mean,na.rm=T)
					if (quitomono==TRUE){
                        #solo quito monomorficos
                        vect2=which(tfreq2$pest!=1 & tfreq2$pest1!=1 & tfreq2$pest!=0 & tfreq2$pest1!=0)
                        if (length(vect2)!=0){
                          tfreq2=tfreq2[vect2,]
                          nummark[[i]]=length(vect2)
                          tfreq2=t(tfreq2[,1:(ncol(tfreq2)-2)])
                          heW[[i]]=heiter(tfreq2)
                        }else{
                             nummark[[i]]=0
                             heW[[i]]="All markers monomorphics"
					    	 outheT[[i]]=which(groups[,1]==i)
                             }
				    }else{
					    vect2=which(tfreq2$pest!=1 & tfreq2$pest1!=1 & tfreq2$pest!=0 & tfreq2$pest1!=0)						
						if (length(vect2)!=0){
						    writemono=cbind(rownames(tfreq2)[-vect2],tfreq2$pest[-vect2])
						    write.csv(writemono,paste0(getwd(),"/Output_MarkMonoGroups/MarkersMonomorphicsForGroup_",selv,"_in level_",levels(pp[,1])[i],".csv"),row.names=F)
						}
                        tfreq2=t(tfreq2[,1:(ncol(tfreq2)-2)])
						nummark[[i]]=nmark
                        heW[[i]]=heiter(tfreq2)  
					}
					
				}
            }
                        
            outheT=unlist(outheT)
			heW=unlist(heW)
            nummark=c(dim(tfreq)[2],unlist(nummark))
			if(length(outheT)!=0){
			   heT=heiter(tfreq[,-outheT]) 							  ## Population He               
			   heWe=heW[-which(is.na(as.numeric(heW))==TRUE)]			   
               heB=heT-mean(as.numeric(heWe),na.rm=TRUE) 	  ## Among individuals whitin groups He 			   			   
			}else{
                heT=HE 							  ## Population He                        
                heB=heT-mean(heW,na.rm=TRUE) 	  ## Among individuals whitin groups He 
			}
            Fst=heB/heT 								  ## Diversity among groups (estadistico F wright el grado de 			
                                                          ## diferencia genetica entre las poblaciones, en funcion de las frecuencias alelicas
            		
            div2=cbind(t(t(c("Diversity among groups","Diversity within group",rep("",(length(heW)-1))))),t(t(c("",levels(as.factor(pp[,1]))))),t(t(c(Fst,heW))))
            div2=cbind(div2,nummark)
            agc.env=as.data.frame(groups[,1])
			names(agc.env)<-c("Pop")
			agc.env$Pop<-as.factor(agc.env$Pop)
			tabamv=forAMOVA(mdrMAT,agc.env)
            if("div2"%in%ls()==TRUE){				
				cat("\n","Population Structure for group: ",selv,"\n","\n",file=out,append=T)			
				write.table(div2, file = out, append = T,quote=F, sep=",",col.names=F,row.names=F)				
				cat("\n","AMOVA for group: ",selv,"\n","\n",file=out,append=T)
				write.table(tabamv, file = out, append = T,quote=F, sep=",",col.names=T,row.names=F)			
            }			
	   }
return(list(div2,tabamv))
}