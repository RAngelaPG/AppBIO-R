## ----setup, include = FALSE---------------------------------------------------
library(bslib)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVWI0cA7lFrlGqqfwCC6diIAEnk7U8BeT11uRWVDKABzOA8Ib29URjhqJRk4RgDPOk5Sdm5SGGoM0kiAGVpsrlQ+AH0ARnYBXE8hMBKUvgjPDG6Wxs9iZPTAloBiAGYCAA4ZegJe3C847jgoVKGg+RD9GQAJFbWGhdi4zlpU+hYM4NCDAGUzuAvGQ8XvRRlWK82bmQAhIg+MTicSgBFItBIvxMpFIJA8YESyVWLSaLSwSRSvVenka2PBpGocAyLUMiJSaRai1xx04aQAbmkMgAzEwQMEQiDsPioaFNIjQnmkJq0zinEgCTwgbEjTxMqAAayJZORx28RHotMYDIAogyyFyIIKACTKmRNKWqk5wUUck0Y1YNbEAX1U3idIgEYCdAF0gA", 
    height = "800", 
    width = "100%", 
    style = "border: 1px solid rgba(0,0,0,0.175); border-radius: .375rem;", 
    allowfullscreen = "", 
    allow = "autoplay", 
    `data-external` = "1"
  ),
  full_screen = TRUE,
  style = "margin: 0 auto; float: none;"
)

## ----eval=FALSE---------------------------------------------------------------
# library(shiny)
# library(shinydashboard)
# library(shinydashboardPlus)
# library(waiter)
# shinyApp(
#   ui = dashboardPage(
#     preloader = list(html = tagList(spin_1(), "Loading ..."), color = "#3c8dbc"),
#     header = dashboardHeader(),
#     sidebar = dashboardSidebar(),
#     body = dashboardBody(
#       actionButton("reload", "Reload")
#     ),
#     title = "Preloader"
#   ),
#   server = function(input, output, session) {
#     # fake reload
#     observeEvent(input$reload, {
#       session$reload()
#     })
#   }
# )

