# options(shiny.error = NA)
library(shiny)
library(shinydashboard)
library(idiogramFISH)
library(rhandsontable)
library(gtools)
library(knitr)
library(rclipboard)
library(clipr)
library(shinyjs)

devel <- FALSE
sep <- ":"
if (Sys.info()["sysname"] == "Windows") {
  sep <- ";"
}

if (system.file(package = "rmarkdown") != "") {
  if (rmarkdown::pandoc_available()) {
    if (rmarkdown::pandoc_version() > 2.11 && devel == FALSE) {
      Sys.setenv(PATH = paste0(Sys.getenv("PATH"), sep, sub(".pandoc$", "", rmarkdown::pandoc_exec())))
      rmarkdown::render("www/README2.Rmd")
    }
  }
}

ui <- tagList(
  tags$style(type = "text/css", "#ARCImessageOut {white-space: pre-wrap; word-break: keep-all;}"),
  tags$style("
  body {
    -moz-transform: scale(0.8, 0.8); /* Moz-browsers */
    zoom: 0.8; /* Other non-webkit browsers */
    zoom: 80%; /* Webkit browsers */
    }
  "),
  tags$style("html,body{background-color: white;}
                .shiny-split-layout {
                  width: 93% !important;
                }
                .container{
                    width: 100%;
                    margin: 0 auto;
                    padding: 0;
                }
               @media screen and (min-width: 700px){
                .container{
                    min-width: 1850px;
                    max-width: 1920px;
                }
               }
             .form-control {
  padding-top: 0px;
  padding-right: 0px;
  padding-bottom: 0px;
  padding-left: 5px;
}
"),
  tags$div(class = "container",
    dashboardPage(title = "idiogramFISH",
      dashboardHeader(title = p("idiogramFISH shiny app")),
      dashboardSidebar(
        uiOutput("mysidebar")
      ),
      dashboardBody(
        rclipboardSetup(),
        useShinyjs(),
        tags$head(
          tags$style(HTML(".fa{font-size: 12px;}")),
          tags$style(HTML(".innerI .fa {font-size: 18px;}")),
          tags$style("#markdown{height:100vh !important;}"),
          tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css"),
          tags$style(
            "
body{
    height: auto;
    margin: auto;
    overflow-x: auto;
}"
          )
        ) # head
        , uiOutput("myTabs")
      ) # dashboardbody
    )
  )
)
