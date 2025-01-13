.libPaths("./R-Portable/App/R-Portable/library")
message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
shiny::runApp("./shiny/",port=8888,launch.browser=TRUE)


