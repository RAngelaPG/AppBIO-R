pks<-c("adegenet","ape","cluster","car","chromoMap","corehunter","data.table","dendextend","DT","dplyr","factoextra","ggplot2","BiocManager","Hmisc","htmlwidgets","naturalsort","pedigreemm","plotly","plyr",
"rJava","RColorBrewer","reshape", "reshape2", "shiny", "shinyBS","shinyFiles","shinyjs","shinyWidgets","shinydashboard","shinydashboardPlus","shinyalert", "stringr", "vcfR", "vegan", "ff")
checkpack<-pks%in%rownames(installed.packages())
if(length(checkpack)!=0) {
  pksinst=pks[which(checkpack==FALSE)]
  install.packages(pksinst,lib=.libPaths()[1])
}
if("ggtree"%in%rownames(installed.packages())==FALSE){
  BiocManager::install("ggtree")
}
install_github("whweve/IntAssoPlot")
#Cambiar este directorio de trabajo "C:/Users/RAPACHECO/OneDrive - CIMMYT/Desktop" por el de cada usuario
shiny::runApp("C:/Users/RAPACHECO/OneDrive - CIMMYT/Desktop/AppBIO-R/shiny",port=8887,launch.browser=TRUE)




