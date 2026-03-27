## ----include = FALSE----------------------------------------------------------
library(bslib)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVVP4BBdOxEACRydr2AvPd3dFyw1ADmcHYQjiGcBIxE1NQAKkSxqG720VgAqgCiuA4h9txwUDJwjIme3gYAEnkFjOwCmcHZnLQF9CzF8l76MgDKTXAt1TB8tABmtHAyicnpuPbEUWic45OpabVZIYoyrG16ygBCRFtB2Y7UaKjUHACMiABMAAwzigAeMwDuTaTciVe3M6S0UjUOCJIRgF5g1QnNb1ELEMgRaj9HYdZQAYRIpER-RqdROAKBIPcYIAIu1SkZ-HAwesYaFCgA3QqJYYmCAEAEkdh8VAmUgzIh83mkAT2ED2AC+DgEYAlAF0gA", 
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
# shinyApp(
#    ui = dashboardPage(
#      options = list(sidebarExpandOnHover = TRUE),
#      header = dashboardHeader(),
#      sidebar = dashboardSidebar(minified = TRUE, collapsed = TRUE),
#      body = dashboardBody(
#       lapply(1:20, box, width = 12, title = "box")
#      ),
#      controlbar = dashboardControlbar(),
#      title = "DashboardPage"
#    ),
#    server = function(input, output) { }
#  )

