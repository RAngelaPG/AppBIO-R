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
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVQAIp-AILp2IixZO0LAXgu7ui5YagBzOAcIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0MUAD2DQ0M4iOhlo2Pj3CwAxGwAZDIBRAsanUlpSajguoTAM0ihSMysTGBgJWZHR3IIAa39GIhMIGS6AOQBVPr6txoB3bNJuLoAWW6aFpc4ZsDkIQMYm2qoQAZkQiORGF0QaZslgiHcGqMnMRTDAIEjkU4HjInl0AGzvZFxTgERi0VDjEilahEPaYrFOCCreidDyzACMAHYAKSbCzAxnMmCsxgAYRaRCh7LABzgFH5gqxwtFAEliCEPLQNQ5CCw4KQALQmVCzJJK5HcGJxaUCsAAEgAzABWXAAJg5AAYMC9HYrCoyxnBaqQfgAVADyYf6FiwgwAaoNLoN-YHQuT-NxSOUDGyLGGsBdhhbRutGP4+OVSKQiDAur0BimA1jLC2iaFUasMSXsY9nh5CT3PHBSeTKbRqbT6UOmSy85yABx8vAC5tCufiyW22aMOAyVNpiwq+LqkhdbUkXUEfVGmQIiBm9uW63zsAct2ez0HtPkEM-ADiEYxmKEYALKGH0gxhqqEZnN+gYZlmOY2vW-RDE+pYsBWEBVjWdYeA2QxDq2yIkaEZEWGR5rNuMkzTDKAAiI5khSVIhDSdK7LMKTUU4nDxAAbnmILHAQbHsHwqAmKQuAWEcpBSaQAgWCAFgAL6OAIYBqQAukAA", 
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
#  shinyApp(
#   ui = dashboardPage(
#     dashboardHeader(),
#     dashboardSidebar(),
#     dashboardBody(
#      box(
#       solidHeader = FALSE,
#       title = "Status summary",
#       background = NULL,
#       width = 4,
#       status = "danger",
#       footer = fluidRow(
#         column(
#           width = 6,
#           descriptionBlock(
#             number = "17%",
#             numberColor = "green",
#             numberIcon = icon("caret-up"),
#             header = "$35,210.43",
#             text = "TOTAL REVENUE",
#             rightBorder = TRUE,
#             marginBottom = FALSE
#           )
#         ),
#         column(
#           width = 6,
#           descriptionBlock(
#             number = "18%",
#             numberColor = "red",
#             numberIcon = icon("caret-down"),
#             header = "1200",
#             text = "GOAL COMPLETION",
#             rightBorder = FALSE,
#             marginBottom = FALSE
#           )
#         )
#       )
#      )
#     ),
#     title = "Description Block"
#   ),
#   server = function(input, output) { }
#  )

