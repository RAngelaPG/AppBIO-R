#Install JDK from java and in the terminal write R CMD javareconf
#Works in R.4.1.2
pks<-c("adegenet","ape","cluster","car","chromoMap","corehunter","data.table","dendextend","DT","dplyr","factoextra","ggplot2","BiocManager","Hmisc","htmlwidgets","naturalsort","pedigreemm","plotly","plyr",
"rJava","RColorBrewer","reshape", "reshape2", "shiny", "shinyBS","shinyFiles","shinyjs","shinyWidgets","shinydashboard","shinydashboardPlus","shinyalert", "stringr", "vcfR", "vegan", "ff","statgenGWAS")
checkpack<-pks%in%rownames(installed.packages())
if(length(checkpack)!=0) {
  pksinst=pks[which(checkpack==FALSE)]
  install.packages(pksinst,lib=.libPaths()[1])
}
if("ggtree"%in%rownames(installed.packages())==FALSE){
  BiocManager::install("ggtree")
  BiocManager::install("rtracklayer")
}
library(devtools)
install_github("whweve/IntAssoPlot")
#Change work directory
shiny::runApp("~/AppBIO-R/shiny",port=8887,launch.browser=TRUE)




