## ----include = FALSE----------------------------------------------------------
library(bslib)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----controlbar-basis-code, eval=FALSE----------------------------------------
# library(shiny)
# library(shinydashboard)
# library(shinydashboardPlus)
# shinyApp(
#   ui = dashboardPage(
#     header = dashboardHeader(),
#     sidebar = dashboardSidebar(),
#     body = dashboardBody(),
#     controlbar = dashboardControlbar(collapsed = FALSE),
#     title = "DashboardPage"
#   ),
#   server = function(input, output) { }
# )

## ----controlbarMenu-code, eval=FALSE------------------------------------------
# menu <- controlbarMenu(
#   id = "controlbarMenu",
#   controlbarItem(
#     "Tab 1",
#     "Welcome to tab 1"
#   ),
#   controlbarItem(
#     "Tab 2",
#     numericInput("num", "Observations:", 200, min = 1, max = 1000, step = 100)
#   )
# )
# 
# shinyApp(
#   ui = dashboardPage(
#     header = dashboardHeader(),
#     sidebar = dashboardSidebar(),
#     body = dashboardBody(),
#     controlbar = dashboardControlbar(
#       skin = "dark",
#       menu
#     ),
#     title = "Right Sidebar"
#   ),
#   server = function(input, output) { }
# )

## ----eval=FALSE---------------------------------------------------------------
# sidebarMenu(
#   id = "sidebarMenu",
#   lapply(1:3, function(i) {
#     menuItem(
#       sprintf("Menu %s", i),
#       tabName = sprintf("menu_%s", i),
#       icon = icon("circle")
#     )
#   })
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSPgQTAAgA8AWhvEyjItXosAshRPsRNja0MjYAvDZChCSkbh7evpG4Ac7RsZ6MAJLkMP4QgYGRACpQ9DYAjInJBWAA6nDUxPA2pETNJeWRyQJJeSmu7ulZcDlVEWDFpQBMlb2B1vCMtAQZEKgmpP6QJjCJYwDy9JxwjABuUKS0JJyIu5MADHe4NjB8YeVPMFAAHm9lD482TjkVC-B6qQKqCwQKT8ACC6FygRMtDeum4imUhigAHM4Ij8tw4FAZMdUfJ0foZAAJIkkxjsbqjTjBODpMl6ZQAZRZ6Xx+WefFoADNaHAQuFClgAKoAUSeo0CxGo1DQR3FNklsp6-MCzJJ6R81j5OuCb0ietZ8WsMx1+RV6GoHDKiAAzE8hSYIAQLiR2LQBDYQAqdVYTEMRrNbflOKhFmQhZtDbYAKScXb++WRqNtegAOVgcDeMbjpATkVDAH1U+nuhEs1GliQ3o2IJsCLRGARqHBIuDszY+7aAL6D-KjxmRxQyVjsinKABCRGnxsCpBK4c4K-ya-o4fYO-zTXC5d8FYqeDGAGFohRSB0wLXg6v19kt-yDwWzWBK9M8E-8ooXxvjqFykN2X6wtQADuUCsJwzSEoCn6oNQRCkAAhDa-aKu49rMvQ4EStKcp1thgQoWhezrGsGyRDItBAsYaG9tqZEWmy4SAdy+osMBtqmseYBcTyLBYWRNieAQADW2JuJ66qRAAxAAHHcqmqWJbF0HSKw0XxUaREQhy7P+2aRLm2z0KSRBCjYRlHKc5yXNCNx-vW2EvHk4SPKZUafD84R-A8mbifyZymIW4QAKwPL5-KjtmCU6klELBuOrE4f0cSMLOGIGNeWW8sGAljC4MQDKJbm2qGowTu+tBgZFYwACLknlRg4j2YBdKxDknKS4Qel6PqtnwNFPEQ1HrE8RycMyJABkGkaTaQNEACT0YxqF3o4NiMBQdJMRsS22rwQLsIwEBEIwORjesa32QIo4jqMik2JyUENQQ3B9OV2XPL4EnyGKdl5OxLBBNkNjhSYcAYDYXjXXARD9Ywr1BLZpCIeDOWhpDwxBPBkxPFjhZlWkENEKgFBpsthzHP1Mr9WQfqrPdONJk8J0mjIPy7UCjAxnQGx3aQa0c748pgBWvbAMAZQALoK8AkwK8Vtl+rzYThJMi3-iYqByOQBV-UVUSFZVSVDsGBtG3AJsU4wSZtqkFVOwkF5HN23og+EwRfM9o7o5edDSXw2KgwhhZHC4ITk27+MwM8UBSXA8GkxJRA-DjdnUxAoz2QzcBM7erPrfH2Wc4G6s2GX90VwagOhIJEw2L+Aa2+ccCLl83GWvSkTCTxaMPqML2RgX9OnMXzMi2zYvzFzNuG13nLaccunrFwafzRAUv2bsMNNaLa3zM9Ty0NiV37SsDVvJqMpB701sQAIYBDgrQA", 
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
# shinyApp(
#   ui = dashboardPage(
#     header = dashboardHeader(),
#     sidebar = dashboardSidebar(
#       minified = TRUE,
#       collapsed = TRUE,
#       sidebarMenu(
#         id = "sidebarMenu",
#         lapply(1:3, function(i) {
#           menuItem(
#             sprintf("Menu %s", i),
#             tabName = sprintf("menu_%s", i),
#             icon = icon("circle")
#           )
#         })
#       )
#     ),
#     body = dashboardBody(
#       tabItems(
#         tabItem(tabName = "menu_1", "Content 1"),
#         tabItem(
#           tabName = "menu_2",
#           box(
#             title = "Always the same plot!",
#             collapsible = TRUE,
#             plotOutput("distPlot"),
#             sidebar = boxSidebar(
#               id = "boxSidebar",
#               background = "#808080",
#               width = "50%",
#               sliderInput(
#                 "obs",
#                 "Number of observations:",
#                 min = 0,
#                 max = 1000,
#                 value = 500
#               )
#             )
#           )
#         )
#       )
#     ),
#     controlbar = dashboardControlbar(
#       id = "controlbar",
#       menu
#     ),
#     title = "DashboardPage"
#   ),
#   server = function(input, output, session) {
#     output$distPlot <- renderPlot({
#       hist(rnorm(input$obs))
#     })
#     # Switch controlbar menu based on sidebar item value. Moreover
#     # if the sidebar menu item is 2, the controlbar opens
#     observeEvent(input$sidebarMenu, {
#       idx <- strsplit(input$sidebarMenu, "_")[[1]][2]
#       if (idx == 2) {
#         updateControlbar("controlbar")
#       }
#       updateControlbarMenu("controlbarMenu", selected = idx)
#     })
# 
#     # Clicking on the second controlbar item makes the box sidebar open
#     observeEvent(input$controlbarMenu, {
#       if (input$controlbarMenu == "Tab 2") updateBoxSidebar("boxSidebar")
#     })
# 
#     observeEvent(input$num, {
#       updateSliderInput(session, "obs", value = input$num)
#     }, ignoreInit = TRUE)
# 
#   }
# )

## ----footer-code, eval=FALSE--------------------------------------------------
# shinyApp(
#    ui = dashboardPage(
#      header = dashboardHeader(),
#      sidebar = dashboardSidebar(),
#      body = dashboardBody(),
#      footer = dashboardFooter(left = "Left content", right = "Right content"),
#      title = "DashboardPage"
#    ),
#    server = function(input, output) { }
#  )

