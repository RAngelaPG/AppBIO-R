.libPaths("./R-Portable/App/R-Portable/library")
Sys.setenv(RSTUDIO_PANDOC=paste0(unlist(strsplit(.libPaths(),"/R-Portable/App/R-Portable/library")),"/pandocExe/"))
rmarkdown::find_pandoc(cache=F,dir=Sys.getenv("RSTUDIO_PANDOC"))
message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB
shiny::runApp("./shiny/",port=8888,launch.browser=TRUE)


