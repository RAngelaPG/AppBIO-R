.libPaths("./R-Portable/App/R-Portable/library")
message('library paths:\n', paste('... ', .libPaths(), sep='', collapse='\n'))
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB
shiny::runApp("./shiny/",port=8888,launch.browser=TRUE)


