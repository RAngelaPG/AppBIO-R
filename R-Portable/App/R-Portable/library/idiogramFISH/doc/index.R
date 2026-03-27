## ----cssjs, results="asis", echo=FALSE, message=FALSE, eval=TRUE--------------
# <!-- pkgdown -->
# <!-- jquery --><script src="js/jquery.min.js" crossorigin="anonymous"></script>

myfile <- "js/jquery.min.js"
if (file.exists(myfile)) {
  cat(paste0('<script src="', myfile, '" crossorigin="anonymous"></script> <!-- # -->'))
}
# <!-- clipboard.js --><script src="js/clipboard.min.js"  crossorigin="anonymous"></script>
myfile <- "js/clipboard.min.js"
if (file.exists(myfile)) {
  cat(paste0('<script src="', myfile, '" crossorigin="anonymous"></script>'))
}
myfile <- "js/hideOutput.js"
if (file.exists(myfile)) {
  cat(paste0('<script src="', myfile, '" crossorigin="anonymous"></script>'))
}
# <!-- Font Awesome icons --><link rel="stylesheet" href="css/all.minMod.css"  crossorigin="anonymous">
myfile <- "css/all.minMod.css"
if (file.exists(myfile)) {
  cat(paste0('<link rel="stylesheet" href="', myfile, '"  crossorigin="anonymous">'))
}
# <!-- Bootstrap --><link rel="stylesheet" href="css/bootstrap.minO.css" crossorigin="anonymous">
myfile <- "css/bootstrap.minO.css"
if (file.exists(myfile)) {
  cat(paste0('<link rel="stylesheet" href="', myfile, '"  crossorigin="anonymous">'))
}
# <!-- # <script src="js/bootstrap.min.js"  crossorigin="anonymous"></script> -->
myfile <- "js/bootstrap.min.js"
if (file.exists(myfile)) {
  cat(paste0('<script src="', myfile, '" crossorigin="anonymous"></script> <!-- # -->'))
}
myfile <- "js/pkgdown2.js"
if (file.exists(myfile)) {
  cat(paste0('<script src="', myfile, '"></script> <!-- # -->'))
}

## ----setup, include=FALSE, eval=T---------------------------------------------

# Create myheader.html

line1 <- '<script src="https://kit.fontawesome.com/af0a13599b.js" crossorigin="anonymous"></script>'
line2 <- '<link rel="shortcut icon" href="../man/figures/logo.png" />'
file <- "myheader.html"

if (Sys.info()["sysname"] == "Windows") {

  # check internet response
  res <- !as.logical(system(paste("ping", "www.google.com")))
  if (res) {
    write(line1, file = file)
    write(line2, file = file, append = TRUE)
  }
} else {
  write(line1, file = file)
  write(line2, file = file, append = TRUE)
}

require(idiogramFISH)

knitr::opts_chunk$set(eval = TRUE)

pasteLinks <- function(link) {
  cat(paste0(
    "
</br>
Jupyter interactive version:
<br>
</br>
<table>
<tr>
<td>
[<img src='../man/figures/colab-badge.svg'>](https://colab.research.google.com/github/fernandoroa/IFjupyter/blob/main/",
    note,
    "){target='_blank'}&emsp;
</td>
<td>
&emsp;&emsp;<img src='../man/figures/GitHub-Mark-120px-plus.png' height=25 width=25>&emsp;[Github](https://github.com/fernandoroa/IFjupyter/blob/main/",
    note,
    "){target='_blank'}
</td>
<td>
&emsp;&emsp;[Raw](https://github.com/fernandoroa/IFjupyter/raw/main/",
    note,
    "){target='_blank'}
</td>
</tr>
</table>
</br>
"
  ))
}

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  # version of manual
  v <- sub("Version: ", "", readLines("../DESCRIPTION")[3])
  # v<-tryCatch(suppressWarnings(rvcheck:::check_github_gitlab("ferroao/idiogramFISH", "gitlab")$latest_version), error=function(e) NA )
  pkg <- "idiogramFISH"
  link <- tryCatch(suppressWarnings(badger::badge_custom("Documentation", paste(pkg, v), "cornflowerblue")), error = function(e) NA)
  if (!is.na(link)) {
    svglink <- gsub("\\[|\\]|!|\\(|\\)", "", link)
    manual_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(svglink)), error = function(e) NA)
    if (!is.na(manual_cont)) {
      manual_contFile <- "../man/figures/manual.svg"
      writeLines(manual_cont, con = manual_contFile)
      manual_contFile <- normalizePath(manual_contFile)
      knitr::include_graphics(manual_contFile)
      # cat(paste0("&nbsp;![''](",knitr::include_graphics(manual_contFile),")" ) )
    }
  }
} # rcurl

## ---- echo=F,  results="asis", eval=T-----------------------------------------
img1_path <- "../man/figures/logo.png"
if (file.exists(img1_path)) {
  cat(paste0("<img src=", img1_path, ">"))
}

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----

if (requireNamespace("RCurl", quietly = TRUE)) {
  cranversion <- "https://www.r-pkg.org/badges/version/idiogramFISH"
  cranversion_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(cranversion)), error = function(e) NA)

  if (!is.na(cranversion_cont)) {
    cranversion_contFile <- "../man/figures/cranversion.svg"
    writeLines(cranversion_cont, con = cranversion_contFile)
    cranversion_contFile <- normalizePath(cranversion_contFile)
    # knitr::include_graphics(cranversion_contFile)
    # cat(paste0("![https://CRAN.R-project.org/package=idiogramFISH](",knitr::include_graphics(cranversion_contFile),")" ) )
    cat(paste0("[![CRAN repo](", knitr::include_graphics(cranversion_contFile), ")](https://CRAN.R-project.org/package=idiogramFISH)"))
  } # cran version

  # crandownloads<-"https://cranlogs.r-pkg.org/badges/grand-total/idiogramFISH?color=orange"
  # crandownloads_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(crandownloads) ), error=function(e) NA )

  doibadge <- "https://zenodo.org/badge/DOI/10.5281/zenodo.3579417.svg"
  doibadge_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(doibadge)), error = function(e) NA)

  if (!is.na(doibadge_cont)) {
    doibadge_contFile <- "../man/figures/doibadge.svg"
    writeLines(doibadge_cont, con = doibadge_contFile)
    doibadge_contFile <- normalizePath(doibadge_contFile)

    cat(paste0("&nbsp;[![10.5281/zenodo.3579417](", knitr::include_graphics(doibadge_contFile), ")](https://doi.org/10.5281/zenodo.3579417)"))
  } # doi
} # rcurl

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  # cran version
  v <- "NEWS" # tryCatch(suppressWarnings(rvcheck::check_cran("idiogramFISH")$latest_version), error=function(e) NA )

  link <- tryCatch(suppressWarnings(badger::badge_custom("gitlab", paste(v), "orange", "?logo=gitlab")), error = function(e) NA)
  if (!is.na(link)) {
    svglink <- gsub("\\[|\\]|!|\\(|\\)", "", link)
    news_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(svglink)), error = function(e) NA)
    if (!is.na(news_cont)) {
      news_cont_contFile <- "../man/figures/NEWS.svg"
      writeLines(news_cont, con = news_cont_contFile)
      cat(paste0("[![NEWS](", knitr::include_graphics(news_cont_contFile), ")](https://gitlab.com/ferroao/idiogramFISH/blob/master/NEWS.md){target='_blank'}"))
    }
  }
} # rcurl

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  archivelink <- tryCatch(suppressWarnings(badger::badge_custom("CRAN", "archive", "gray")), error = function(e) NA)
  if (!is.na(archivelink)) {
    svgnewdownlink <- gsub("\\[|\\]|!|\\(|\\)", "", archivelink)
    archive_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(svgnewdownlink)), error = function(e) NA)
  }
  if (!is.na(archive_cont)) {
    archive_contFile <- "../man/figures/archive.svg"
    writeLines(archive_cont, con = archive_contFile)
    cat(paste0("[![archive](", knitr::include_graphics(archive_contFile), ")](https://cran.r-project.org/src/contrib/Archive/idiogramFISH/){target='_blank'}"))
  }
} # rcurl

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  newdownlink <- tryCatch(suppressWarnings(badger::badge_custom("CRAN", "downloads", "green")), error = function(e) NA)
  if (!is.na(newdownlink)) {
    svgnewdownlink <- gsub("\\[|\\]|!|\\(|\\)", "", newdownlink)
    realdo_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(svgnewdownlink)), error = function(e) NA)
  }
  if (!is.na(realdo_cont)) {
    realdo_contFile <- "../man/figures/realdownloads.svg"
    writeLines(realdo_cont, con = realdo_contFile)
    # message                                              # link
    cat(paste0("[![downloads](", knitr::include_graphics(realdo_contFile), ")](https://ferroao.gitlab.io/idiogramfishhelppages/downloads.png){target='_blank'}"))
  }
} # rcurl

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  # version of manual
  # v<-tryCatch(suppressWarnings(rvcheck:::check_github_gitlab("ferroao/idiogramFISH", "gitlab")$latest_version), error=function(e) NA )
  # cran version
  v <- tryCatch(suppressWarnings(rvcheck::check_cran("idiogramFISH")$latest_version), error = function(e) NA)

  pkg <- "idiogramFISH"
  link <- tryCatch(suppressWarnings(badger::badge_custom("bookdown", paste(pkg, v), "orange")), error = function(e) NA)
  if (!is.na(link)) {
    svglink <- gsub("\\[|\\]|!|\\(|\\)", "", link)
    manual_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(svglink)), error = function(e) NA)
    if (!is.na(manual_cont)) {
      bookdown_contFile <- "../man/figures/cranmanualbookdown.svg"
      writeLines(manual_cont, con = bookdown_contFile)
      #    manual_contFile <- normalizePath(manual_contFile)
      # knitr::include_graphics(manual_contFile)
      cat(paste0("[![https://ferroao.gitlab.io/manualidiogramfish](", knitr::include_graphics(bookdown_contFile), ")](https://ferroao.gitlab.io/manualidiogramfish/)"))
    }
  }
} # rcurl

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  # version of manual
  # v<-sub("Version: ","",readLines("DESCRIPTION")[3])
  v <- tryCatch(suppressWarnings(rvcheck:::check_github_gitlab("ferroao/idiogramFISH", "gitlab")$latest_version), error = function(e) NA)
  # cran version
  # v<-tryCatch(suppressWarnings(rvcheck::check_cran("idiogramFISH")$latest_version), error=function(e) NA )

  pkg <- "idiogramFISH"
  pkglink <- tryCatch(suppressWarnings(badger::badge_custom("pkgdown", paste(pkg, v), "cornflowerblue")), error = function(e) NA)
  if (!is.na(pkglink)) {
    pkgsvglink <- gsub("\\[|\\]|!|\\(|\\)", "", pkglink)
    develpkgmanual_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(pkgsvglink)), error = function(e) NA)
    if (!is.na(develpkgmanual_cont)) {
      pkgdownmanual_contFile <- "../man/figures/pkgdownver.svg"
      writeLines(develpkgmanual_cont, con = pkgdownmanual_contFile)
      cat(paste0("[![https://ferroao.gitlab.io/idiogramFISH](", knitr::include_graphics(pkgdownmanual_contFile), " )](https://ferroao.gitlab.io/idiogramFISH)"))
    }
  }
} # rcurl

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis"----
if (requireNamespace("RCurl", quietly = TRUE)) {
  # version of manual
  # v<-sub("Version: ","",readLines("DESCRIPTION")[3])
  v <- tryCatch(suppressWarnings(rvcheck:::check_github_gitlab("ferroao/idiogramFISH", "gitlab")$latest_version), error = function(e) NA)
  # cran version
  # v<-tryCatch(suppressWarnings(rvcheck::check_cran("idiogramFISH")$latest_version), error=function(e) NA )

  pkg <- "idiogramFISH"
  link <- tryCatch(suppressWarnings(badger::badge_custom("vignettes", paste(pkg, v), "cornflowerblue")), error = function(e) NA)
  if (!is.na(link)) {
    vignettelink <- gsub("\\[|\\]|!|\\(|\\)", "", link)
    vignettemanual_cont <- tryCatch(suppressWarnings(RCurl::getURLContent(vignettelink)), error = function(e) NA)
    if (!is.na(vignettemanual_cont)) {
      vignettemanual_contFile <- "../man/figures/develmanualvignette.svg"
      writeLines(vignettemanual_cont, con = vignettemanual_contFile)
      #    manual_contFile <- normalizePath(manual_contFile)
      # knitr::include_graphics(manual_contFile)
      cat(paste0("[![https://ferroao.gitlab.io/idiogramfishhelppages](", knitr::include_graphics(vignettemanual_contFile), ")](https://ferroao.gitlab.io/idiogramfishhelppages)"))
    }
  }
} # rcurl

## ----citation, echo=FALSE, comment=NA,results='asis'--------------------------
print(citation("idiogramFISH"), bibtex = FALSE)

## ---- echo=TRUE, results=FALSE, eval=FALSE------------------------------------
#  sink("idiogramFISH.bib")
#  toBibtex(citation("idiogramFISH"))
#  sink()

## ---- echo=F,  results="asis", eval=T-----------------------------------------
img1_path <- "../man/figures/kofi1.png"
img1_path <- normalizePath(img1_path)

if (file.exists(img1_path)) {
  cat(paste0("[Fernando Roa](https://ferroao.gitlab.io/curriculumpu/){target='_blank'}&nbsp;&nbsp;<a href='https://ko-fi.com/X7X71PZZG' target='_blank'><img src=", img1_path, " width=\"10%\">"))
} else {
  cat(paste0("[Fernando Roa](https://ferroao.gitlab.io/curriculumpu/){target='_blank'}"))
}

## ---- echo=F, message=FALSE, fig.show = "hold", fig.align = "default", results="asis", eval=FALSE----
#  cat("<a href='https://ko-fi.com/X7X71PZZG' target='_blank'><img height='12' style='border:0px;height:12px;' src='../man/figures/kofi1.png' border='0' alt='Buy Me a Coffee at ko-fi.com' /></a>")

## ----include=FALSE,eval=FALSE-------------------------------------------------
#  # automatically create a bib database for R packages, this is currently not used by vignette refs/packages2.bib
#  knitr::write_bib(c(
#    .packages(), "bookdown", "knitr", "rmarkdown", "devtools", "pkgdown", "crayon", "ggtree", "ggplot2", "ggpubr", "phytools", "plyr", "dplyr", "tidyr", "rentrez"
#  ), "refs/packages2.bib")

## ---- echo=F------------------------------------------------------------------
chapterFile_jupyter <- "../parts/01-jupyter.Rmd"
if (file.exists(chapterFile_jupyter)) {
  childExists_jupyter <- TRUE
  child_docs_jupyter <- chapterFile_jupyter
} else {
  childExists_jupyter <- FALSE
  child_docs_jupyter <- ""
}

## ---- echo=F------------------------------------------------------------------
chapterFile4 <- "../chaptersBLOCK/12-human.Rmd"
if (file.exists(chapterFile4)) {
  renderLink <- FALSE
} else {
  renderLink <- TRUE
}

## ---- results="asis", echo=FALSE----------------------------------------------
cat("Shiny App in the cloud availability:  \n")
cat("[shinyapps.io](https://ferapps.shinyapps.io/iboard/){target='_blank'}")

## ---- echo=F,  results="asis", eval=T-----------------------------------------
jupyter_path <- "../man/figures/colab2.jpg"
jupyter_path <- normalizePath(jupyter_path)
if (file.exists(jupyter_path)) {
  cat(paste0("![](", jupyter_path, ")")) # works
}

## ---- echo=F------------------------------------------------------------------
chapterFile_first <- "../parts/02-first.Rmd"
if (file.exists(chapterFile_first)) {
  childExists_first <- TRUE
  child_docs_first <- chapterFile_first
} else {
  childExists_first <- FALSE
  child_docs_first <- ""
}

## ---- echo=F------------------------------------------------------------------
chapterFile4 <- "../chaptersBLOCK/12-human.Rmd"
if (file.exists(chapterFile4)) {
  renderLink <- FALSE
} else {
  renderLink <- TRUE
}

## ---- results="asis", echo=FALSE----------------------------------------------
cat("Shiny App in the cloud availability:  \n")
cat("[shinyapps.io](https://ferapps.shinyapps.io/iboard/){target='_blank'}")

## ---- echo=F,  results="asis", eval=F, out.width='100%'-----------------------
#  shiny_path <- "../man/figures/shiny.jpg"
#  if (file.exists(shiny_path)) {
#    cat(paste0("<img src=", shiny_path, ">"))
#  }

## ---- echo=F,  results="asis", eval=T-----------------------------------------
shiny_path <- "../man/figures/shiny.jpg"
shiny_path <- normalizePath(shiny_path)
if (file.exists(shiny_path)) {
  cat(paste0("![](", shiny_path, ")")) # works
}

## ---- echo=F------------------------------------------------------------------
chapterFile_minimal <- "../chaptersBLOCK/03-minimal.Rmd"
if (file.exists(chapterFile_minimal)) {
  childExists_minimal <- TRUE
  child_docs_minimal <- chapterFile_minimal
} else {
  childExists_minimal <- FALSE
  child_docs_minimal <- ""
}

## ---- results="asis", eval=!childExists_minimal, echo=FALSE-------------------
cat("# Minimal Examples")

## ---- results="asis", eval=!childExists_minimal, echo=FALSE-------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#minimal-examples)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists_minimal, comment=NA----
note <- "03-minimal.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile_plotting <- "../chaptersBLOCK/04-plotting.Rmd"
if (file.exists(chapterFile_plotting)) {
  childExists_plotting <- TRUE
  child_docs_plotting <- chapterFile_plotting
} else {
  childExists_plotting <- FALSE
  child_docs_plotting <- ""
}

## ---- results="asis", eval=!childExists_plotting, echo=FALSE------------------
cat("# Plotting chromosomes")

## ---- results="asis", eval=!childExists_plotting, echo=FALSE------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#plotting-chromosomes)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists_plotting, comment=NA----
note <- "04-plotting.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile_multiple <- "../chaptersBLOCK/05-multiple.Rmd"
if (file.exists(chapterFile_multiple)) {
  childExists_multiple <- TRUE
  child_docs_multiple <- chapterFile_multiple
} else {
  childExists_multiple <- FALSE
  child_docs_multiple <- ""
}

## ---- results="asis", eval=!childExists_multiple, echo=FALSE------------------
cat("# Several OTUs")

## ---- results="asis", eval=!childExists_multiple, echo=FALSE------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#multiple-otus)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists_multiple, comment=NA----
note <- "05-multiple.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile_units <- "../chaptersBLOCK/06-units.Rmd"
if (file.exists(chapterFile_units)) {
  childExists_units <- TRUE
  child_docs_units <- chapterFile_units
} else {
  childExists_units <- FALSE
  child_docs_units <- ""
}

## ---- results="asis", eval=!childExists_units, echo=FALSE---------------------
cat("# Changing Units")

## ---- results="asis", eval=!childExists_units, echo=FALSE---------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#changing-units)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists_units, comment=NA----
note <- "06-units.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile0 <- "../chaptersBLOCK/07-gish.Rmd"
if (file.exists(chapterFile0)) {
  childExists0 <- TRUE
  child_docs0 <- chapterFile0
} else {
  childExists0 <- FALSE
  child_docs0 <- ""
}

## ---- results="asis", eval=!childExists0, echo=FALSE--------------------------
cat("# GISH")

## ---- results="asis", eval=!childExists0, echo=FALSE--------------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#gish)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists0, comment=NA----
note <- "07-gish.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile1 <- "../chaptersBLOCK/08-groups.Rmd"
if (file.exists(chapterFile1)) {
  childExists1 <- TRUE
  child_docs1 <- chapterFile1
} else {
  childExists1 <- FALSE
  child_docs1 <- ""
}

## ---- results="asis", eval=!childExists1, echo=FALSE--------------------------
cat("# Using groups")

## ---- results="asis", eval=!childExists1, echo=FALSE--------------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#groups)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists1, comment=NA----
note <- "08-groups.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile2 <- "../chaptersBLOCK/09-circular.Rmd"
if (file.exists(chapterFile2)) {
  childExists2 <- TRUE
  child_docs2 <- chapterFile2
} else {
  childExists2 <- FALSE
  child_docs2 <- ""
}

## ---- results="asis", eval=!childExists2, echo=FALSE--------------------------
cat("# Circular plots")

## ---- results="asis", eval=!childExists2, echo=FALSE--------------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#circular-plots)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists2, comment=NA----
note <- "09-circular.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile3 <- "../chaptersBLOCK/10-phylogeny.Rmd"
if (file.exists(chapterFile3)) {
  childExists3 <- TRUE
  child_docs3 <- chapterFile3
} else {
  childExists3 <- FALSE
  child_docs3 <- ""
}

## ---- results="asis", eval=!childExists3, echo=FALSE--------------------------
cat("# Plotting alongside phylogeny")

## ---- results="asis", eval=!childExists3, echo=FALSE--------------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#plotting-alongside-phylogeny)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists3, comment=NA----
note <- "10-phylogeny.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile <- "../chaptersBLOCK/11-citrushelp.Rmd"
if (file.exists(chapterFile)) {
  childExists <- TRUE
  child_docs <- chapterFile
} else {
  childExists <- FALSE
  child_docs <- ""
}

## ---- results="asis", eval=!childExists, echo=FALSE---------------------------
cat("# *Citrus*")

## ---- results="asis", eval=!childExists, echo=FALSE---------------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#citrus)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists, comment=NA----
note <- "11-citrushelp.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile4 <- "../chaptersBLOCK/12-human.Rmd"
if (file.exists(chapterFile4)) {
  childExists4 <- TRUE
  child_docs4 <- chapterFile4
} else {
  childExists4 <- FALSE
  child_docs4 <- ""
}

## ---- results="asis", eval=!childExists4, echo=FALSE--------------------------
cat("# Human karyotype")

## ---- results="asis", eval=!childExists4, echo=FALSE--------------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#human-karyotype)")

## ---- echo=FALSE, results="asis", message=FALSE, eval=!childExists4, comment=NA----
note <- "12-human.ipynb"
pasteLinks(note)

## ---- echo=F------------------------------------------------------------------
chapterFile_param <- "../chaptersBLOCK/13-functions.Rmd"
if (file.exists(chapterFile_param)) {
  childExists_param <- TRUE
  child_docs_param <- chapterFile_param
} else {
  childExists_param <- FALSE
  child_docs_param <- ""
}

## ---- results="asis", eval=!childExists_param, echo=FALSE---------------------
cat("# Functions")

## ---- results="asis", eval=!childExists_param, echo=FALSE---------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#functions)")

## ---- echo=F------------------------------------------------------------------
chapterFile_data <- "../chaptersBLOCK/14-datasets.Rmd"
if (file.exists(chapterFile_data)) {
  childExists_data <- TRUE
  child_docs_data <- chapterFile_data
} else {
  childExists_data <- FALSE
  child_docs_data <- ""
}

## ---- results="asis", eval=!childExists_data, echo=FALSE----------------------
cat("# Datasets")

## ---- results="asis", eval=!childExists_data, echo=FALSE----------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#datasets)")

## ---- echo=F------------------------------------------------------------------
chapterFile_news <- "../chaptersBLOCK/15-news.Rmd"
if (file.exists(chapterFile_news)) {
  childExists_news <- TRUE
  child_docs_news <- chapterFile_news
} else {
  childExists_news <- FALSE
  child_docs_news <- ""
}

## ---- results="asis", eval=!childExists_news, echo=FALSE----------------------
cat("# News {-}")

## ---- results="asis", eval=!childExists_news, echo=FALSE----------------------
cat("[https://ferroao.gitlab.io/idiogramfishhelppages](https://ferroao.gitlab.io/idiogramfishhelppages/#news)")

