export.xml<-function(filearch,line1,trees,group=0,colorCodes=NULL,groupCodes=NULL,colorCodesB=NULL,groupCodesB=NULL){  
  j=0  
  k=0
  lab1=trees$edge[which(trees$edge[,2]<=length(trees$tip.label)),2]
  write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",filearch)
  write("<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">",filearch,append=TRUE)
  write("<phylogeny rooted=\"false\" rerootable=\"true\">",filearch,append=TRUE)
  write("<clade>",filearch,append=TRUE)
  
  for(i in 1:length(line1)){
    
    if(line1[i]=="(" && line1[i+1]=="("){
      j=j+1
      write("<clade>",filearch,append=TRUE)
      write(paste("<branch_length>",trees$edge.length[j],"</branch_length>",sep=""),filearch,append=TRUE)
    }
    if(line1[i]=="," && line1[i+1]=="("){
      j=j+1
      write("<clade>",filearch,append=TRUE)
      write(paste("<branch_length>",trees$edge.length[j],"</branch_length>",sep=""),filearch,append=TRUE)
    }
    if(line1[i]=="g" && line1[i-1]==","){
      j=j+1
      k=k+1
      write("<clade>",filearch,append=TRUE)
      write(paste("<name>",trees$tip.label[lab1[k]],"</name>",sep=""),filearch,append=TRUE)
      write(paste("<branch_length>",trees$edge.length[j],"</branch_length>",sep=""),filearch,append=TRUE)
	  if(group==1){
	    if(k==1){
		  write("<color>",filearch,append=TRUE)
		  write("<red>255</red>",filearch,append=TRUE)
		  write("<green>255</green>",filearch,append=TRUE)
		  write("<blue>255</blue>",filearch,append=TRUE)
          write("</color>",filearch,append=TRUE)		  
		}
		write(paste("<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">",colorCodes[groupCodes[lab1[k]]],"</property>",sep=""),filearch,append=TRUE)
	  }  
	  if(group==2) {
	     write("<color>",filearch,append=TRUE)
		 write(paste("<red>",colorCodesB[1,groupCodesB[lab1[k]]],"</red>",sep=""),filearch,append=TRUE)
		 write(paste("<green>",colorCodesB[2,groupCodesB[lab1[k]]],"</green>",sep=""),filearch,append=TRUE)
		 write(paste("<blue>",colorCodesB[3,groupCodesB[lab1[k]]],"</blue>",sep=""),filearch,append=TRUE)		 
		 write("</color>",filearch,append=TRUE)
	     write(paste("<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">",colorCodes[groupCodes[lab1[k]]],"</property>",sep=""),filearch,append=TRUE)
	  }      
      write("</clade>",filearch,append=TRUE)
    }
    if(line1[i]=="(" && line1[i+1]=="g"){
      j=j+1
      k=k+1
      write("<clade>",filearch,append=TRUE)
      write(paste("<name>",trees$tip.label[lab1[k]],"</name>",sep=""),filearch,append=TRUE)
      write(paste("<branch_length>",trees$edge.length[j],"</branch_length>",sep=""),filearch,append=TRUE)
	  if(group==1){
	    if(k==1){
		  write("<color>",filearch,append=TRUE)
		  write("<red>255</red>",filearch,append=TRUE)
		  write("<green>255</green>",filearch,append=TRUE)
		  write("<blue>255</blue>",filearch,append=TRUE)
          write("</color>",filearch,append=TRUE)		  
		}
		write(paste("<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">",colorCodes[groupCodes[lab1[k]]],"</property>",sep=""),filearch,append=TRUE)
	  }	  
	  if(group==2) {
	     write("<color>",filearch,append=TRUE)
		 write(paste("<red>",colorCodesB[1,groupCodesB[lab1[k]]],"</red>",sep=""),filearch,append=TRUE)
		 write(paste("<green>",colorCodesB[2,groupCodesB[lab1[k]]],"</green>",sep=""),filearch,append=TRUE)
		 write(paste("<blue>",colorCodesB[3,groupCodesB[lab1[k]]],"</blue>",sep=""),filearch,append=TRUE)		 
		 write("</color>",filearch,append=TRUE)
	     write(paste("<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">",colorCodes[groupCodes[lab1[k]]],"</property>",sep=""),filearch,append=TRUE)
	  }  
      write("</clade>",filearch,append=TRUE)
    }
    if(line1[i]==")"){
        write("</clade>",filearch,append=TRUE)
      }
    
  }
  write("</phylogeny>",filearch,append=TRUE)
  write("</phyloxml>",filearch,append=TRUE)  
}