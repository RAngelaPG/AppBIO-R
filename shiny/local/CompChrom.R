CompChrom<-function(posit,posit1,filename1,filename2){


posit$POS=as.numeric(posit$POS)
posit$CHROM=as.character(posit$CHROM)
tabamv=data.frame(table(as.factor(posit$CHROM)))
colnames(tabamv)=c("Chromosome","NMarks")
chromi=as.character(tabamv[which(tabamv$NMarks>=(sum(tabamv$NMarks)*0.05)),1])
chromi11=chromi[grep("chr",chromi)]
if(length(chromi11)==0){chromi=chromi}else{chromi=chromi11}
if(length(which(posit$CHROM%in%chromi==F))!=0){posit$CHROM[which(posit$CHROM%in%chromi==F)]=NA}		
blast=which(is.na(posit[,1]))

posit1$POS=as.numeric(posit1$POS)
posit1$CHROM=as.character(posit1$CHROM)
tabamv1=data.frame(table(as.factor(posit1$CHROM)))
colnames(tabamv1)=c("Chromosome","NMarks")
chromi1=as.character(tabamv1[which(tabamv1$NMarks>=(sum(tabamv1$NMarks)*0.05)),1])
#chromi1=chromi1[grep("chr",chromi1)]
chromi12=chromi1[grep("chr",chromi)]
if(length(chromi12)==0){chromi1=chromi1}else{chromi1=chromi12}
if(length(which(posit1$CHROM%in%chromi1==F))!=0){posit1$CHROM[which(posit1$CHROM%in%chromi1==F)]=NA}		
blast1=which(is.na(posit1[,1]))

if(length(blast)<length(blast1)){
  positboth=merge(posit,posit1,"ID")
  tit1=filename1
  tit2=filename2
}else{
  positboth=merge(posit1,posit,"ID")
  tit1=filename2
  tit2=filename1
}

c1=str_locate(positboth$CHROM.x, pattern = "chr")[,1]
if(all(is.na(c1))==T){
  c1=positboth$CHROM.x
}else{
  c1=substring(positboth$CHROM.x, first = (c1-2)) 
}

c2=str_locate(positboth$CHROM.y, pattern = "chr")[,1]
if(all(is.na(c2))==T){
  c2=positboth$CHROM.y
}else{
  c2=substring(positboth$CHROM.y, first = (c2-2)) 
}

c11=str_sort(unique(c1), numeric = TRUE)
c21=str_sort(unique(c2), numeric = TRUE)
c11=c11[!is.na(c11)]
c21=c21[!is.na(c21)]

nump=c(paste0("0",1:9),10:1000)
nump=nump[1:length(c11)]
nump=paste0("Chr",nump)

for ( i in 1:length(nump)){
  c1[which(c1==c11[i])]=nump[i]
  c2[which(c2==c21[i])]=nump[i]
}
positboth$CHROM.x=c1
positboth$CHROM.y=c2

positboth$group="both"
positboth$match="match"
positboth$CHROM=positboth$CHROM.x
positboth$POS=positboth$POS.x
for(i in 1:dim(positboth)[1]){
  if(!is.na(positboth$CHROM.x[i]) & is.na(positboth$CHROM.y[i])){positboth$group[i]="g1";positboth$CHROM[i]=positboth$CHROM.x[i];positboth$POS[i]=positboth$POS.x[i];positboth$match[i]="unknown"}
  if(!is.na(positboth$CHROM.y[i]) & is.na(positboth$CHROM.x[i])){positboth$group[i]="g2";positboth$CHROM[i]=positboth$CHROM.y[i];positboth$POS[i]=positboth$POS.y[i];positboth$match[i]="unknown"}
  if(is.na(positboth$CHROM.y[i]) & is.na(positboth$CHROM.x[i])){positboth$group[i]="delete"}
  if(!is.na(positboth$CHROM.y[i]) & !is.na(positboth$CHROM.x[i])){if(positboth$CHROM.y[i]!=positboth$CHROM.x[i]){positboth$match[i]="NOmatch"}}
}

if(length(which(positboth$group=="delete"))!=0){positboth=positboth[-which(positboth$group=="delete"),]}


positjoint=positboth[,c("ID","CHROM","POS","group")]
names(positjoint)=c("ID","CHROM","POS","group")
positjoint$CHROM=as.character(positjoint$CHROM)
positjoint$POS=as.numeric(positjoint$POS)
positjoint$group=as.character(positjoint$group)


chromosome_file=aggregate(POS ~ CHROM, data = positjoint, max)
chromosome_file=data.frame(cbind(CHROM=chromosome_file$CHROM,STAR=1,END=chromosome_file$POS))
chromosome_file$STAR=as.numeric(chromosome_file$STAR)
chromosome_file$END=as.numeric(chromosome_file$END)

matches.x=aggregate(positboth$match, by=list(positboth$CHROM.x,positboth$match), FUN=length)
matches.x=cast(matches.x,Group.1~Group.2,sum)
names(matches.x)[1]="Var1"
revmatch.x=which(c("Var1","match","NOmatch","unknown")%in%names(matches.x)==T)
tmp=data.frame(matrix(0,dim(matches.x)[1],4))
if(length(revmatch.x)!=0){
  for(mt in 1:length(revmatch.x)){
    tmp[,revmatch.x[mt]]=matches.x[,mt]
  }
  matches.x=tmp
  names(matches.x)=c("Var1","match","NOmatch","unknown")
}

matches.y=aggregate(positboth$match, by=list(positboth$CHROM.y,positboth$match), FUN=length)
matches.y=cast(matches.y,Group.1~Group.2,sum)
names(matches.y)[1]="Var1"
revmatch.y=which(c("Var1","match","NOmatch","unknown")%in%names(matches.y)==T)
tmpy=data.frame(matrix(0,dim(matches.y)[1],4))
if(length(revmatch.y)!=0){
  for(mt in 1:length(revmatch.y)){
    tmp[,revmatch.y[mt]]=matches.y[,mt]
  }
  matches.y=tmp
  names(matches.y)=c("Var1","match","NOmatch","unknown")
}

matches=merge(matches.x,matches.y,"Var1")

final=merge(data.frame(table(as.factor(positboth$CHROM.x))),
            data.frame(table(as.factor(positboth$CHROM.y))),"Var1")
final=merge(final,matches,"Var1")
names(final)=c("CHROM",filename1,filename2,paste0("Match_",filename1),paste0("NOmatch_",filename1),paste0("Unknown_",filename1),paste0("Match_",filename2),paste0("NOmatch_",filename2),paste0("Unknown_",filename2))

out=paste0("CompareBlast_",filename1,"_withBlast_",filename2,".csv")
cat("\n","Markers in blast comparison for each chromosome: ","\n","\n",file=out,append=T)
write.table(final, file = out, append = T,quote=F, sep=",",col.names=T,row.names=F)
system2('open',args=out,wait=F)


gg1=which(positboth$group=="g1")
gg2=which(positboth$group=="g2")
ggb=which(positboth$group=="both")

annotation_file1=data.frame(cbind(ID=positboth[c(gg1,ggb),"ID"],CHROM=positboth[c(gg1,ggb),"CHROM.x"],STAR=positboth[c(gg1,ggb),"POS.x"],END=positboth[c(gg1,ggb),"POS.x"]+10))
annotation_file1$STAR=as.numeric(annotation_file1$STAR)
annotation_file1$END=as.numeric(annotation_file1$END)

annotation_file2=data.frame(cbind(ID=positboth[c(gg2,ggb),"ID"],CHROM=positboth[c(gg2,ggb),"CHROM.y"],STAR=positboth[c(gg2,ggb),"POS.y"],END=positboth[c(gg2,ggb),"POS.y"]+10))
annotation_file2$STAR=as.numeric(annotation_file2$STAR)
annotation_file2$END=as.numeric(annotation_file2$END)

return(list(chromosome_file,annotation_file1,annotation_file2))

}