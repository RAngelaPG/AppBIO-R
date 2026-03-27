## ---- echo=F, warning=FALSE, error=FALSE, comment=NA--------------------------
if(Sys.info()['sysname']=="Darwin") {
      system('echo "---" > index.Rmd')
      system('echo "title: \'idiogramFISH\'" >> index.Rmd')
      system('echo "author: \'Fernando Roa\'" >> index.Rmd')
      system('echo "date: \'07 04 2023\'" >> index.Rmd')
      system('echo "output:" >> index.Rmd')
      system('echo "  html_document" >> index.Rmd')
      system('echo "vignette: >" >> index.Rmd')
      system('echo "  %\\VignetteIndexEntry{idiogramFISH}" >> index.Rmd')
      system('echo "  %\\VignetteEngine{knitr::rmarkdown}" >> index.Rmd')
      system('echo "  %\\VignetteEncoding{UTF-8}" >> index.Rmd')
      system('echo "---" >> index.Rmd')
      system('echo ""    >> index.Rmd')
      system('echo "visit https://ferroao.gitlab.io/idiogramfishhelppages" >> index.Rmd')
}
if(length( rmarkdown::pandoc_version() < "2" ) > 0) {
  
  if(rmarkdown::pandoc_version() < "2") {
    
    message(crayon::red("\nMissing pandoc version > 2. Vignette may fail because it uses lua filter for multiple bibliographies
                        \nMore info:
                        \nhttps://stat.ethz.ch/pipermail/r-package-devel/2019q2/004127.html
                        \nhttps://stat.ethz.ch/pipermail/r-package-devel/2020q1/004814.html
                        \nhttps://cran.r-project.org/web/packages/rmarkdown/vignettes/lua-filters.html
                        \nhttps://rmarkdown.rstudio.com/docs/reference/pandoc_args.html
                        \nhttps://pandoc.org/lua-filters.html
                        \nhttps://github.com/jgm/pandoc/releases")
            )
  }
}
if(Sys.info()['sysname']=="Windows") {
    message("pandoc > 2 not available, see online vignettes")
    shell("del index.Rmd")
    
    #
    #    create new index.Rmd
    #
      shell("@echo off")
      shell('@echo --- > index.Rmd')
      shell('@echo title: "idiogramFISH" >> index.Rmd')
      shell('@echo author: "Fernando Roa" >> index.Rmd')
      shell('@echo date: "07 04 2023" >> index.Rmd')
      shell('@echo output: >> index.Rmd')
      shell('@echo   html_document >> index.Rmd')
      shell('@echo vignette: ^> >> index.Rmd')
      shell("@echo   %\\VignetteIndexEntry{idiogramFISH} >> index.Rmd")
      shell('@echo   %\\VignetteEngine{knitr::rmarkdown} >> index.Rmd')
      shell('@echo   %\\VignetteEncoding{UTF-8} >> index.Rmd')
      shell('@echo --- >> index.Rmd')
      shell('@echo(     >> index.Rmd')
      shell('@echo visit https://ferroao.gitlab.io/idiogramfishhelppages >> index.Rmd')
   
} # if windows

