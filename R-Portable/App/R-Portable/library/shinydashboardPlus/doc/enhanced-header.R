## ----setup, include = FALSE---------------------------------------------------
library(bslib)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----leftUi-code, eval=FALSE--------------------------------------------------
# library(shiny)
# library(shinyWidgets)
# library(shinydashboard)
# library(shinydashboardPlus)
#  shinyApp(
#    ui = dashboardPage(
#      header = dashboardHeader(
#        leftUi = tagList(
#          dropdownButton(
#            label = "Controls",
#            icon = icon("sliders-h"),
#            status = "primary",
#            circle = FALSE,
#            sliderInput(
#              inputId = "n",
#              label = "Number of observations",
#              min = 10, max = 100, value = 30
#            ),
#            prettyToggle(
#              inputId = "na",
#              label_on = "NAs kept",
#              label_off = "NAs removed",
#              icon_on = icon("check"),
#              icon_off = icon("trash")
#            )
#          ),
#          dropdownMenu(
#            type = "messages",
#            badgeStatus = "success",
#            messageItem(from = "Support Team", message = "This is the content of a message.", time = "5 mins"),
#            messageItem(from = "Support Team", message = "This is the content of another message.", time = "2 hours"),
#            messageItem(from = "New User", message = "Can I get some help?", time = "Today")
#          )
#        )
#      ),
#      sidebar = dashboardSidebar(),
#      body = dashboardBody(
#        setShadow(class = "dropdown-menu")
#      ),
#      title = "DashboardPage"
#    ),
#    server = function(input, output) { }
#  )

## ----dropdownBlock-code, eval=FALSE-------------------------------------------
# shinyApp(
#   ui = dashboardPage(
#     header = dashboardHeader(
#       leftUi = tagList(
#         dropdownBlock(
#           id = "mydropdown",
#           title = "Dropdown 1",
#           icon = "sliders-h",
#           sliderInput(
#             inputId = "n",
#             label = "Number of observations",
#             min = 10, max = 100, value = 30
#           ),
#           prettyToggle(
#             inputId = "na",
#             label_on = "NAs kept",
#             label_off = "NAs removed",
#             icon_on = icon("check"),
#             icon_off = icon("trash")
#           )
#         ),
#         dropdownBlock(
#           id = "mydropdown2",
#           title = "Dropdown 2",
#           icon = "sliders-h",
#           prettySwitch(
#             inputId = "switch4",
#             label = "Fill switch with status:",
#             fill = TRUE,
#             status = "primary"
#           ),
#           prettyCheckboxGroup(
#             inputId = "checkgroup2",
#             label = "Click me!",
#             thick = TRUE,
#             choices = c("Click me !", "Me !", "Or me !"),
#             animation = "pulse",
#             status = "info"
#           )
#         )
#       )
#     ),
#     sidebar = dashboardSidebar(),
#     body = dashboardBody(
#       setShadow(class = "dropdown-menu")
#     ),
#     title = "DashboardPage"
#   ),
#   server = function(input, output) { }
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJyduOCgZOEZ3T3lvfRkACQioxnYzaIB5E1JUXPswTMYhMAEBXEdQzloo+hZYrx8DAGVauHr0iqqnRRlWRvjmmQAhIn72bpDQ0lpSajhY0oBVTmjPRiJUGSIAdwhSqqmnNcYAN3WPADMTCAJZknY+fNJcByJcl4EHEB6PvNyABJig4ADwAWgcjAoaVW0XYv2mTiaiTh6R6oQc0HgSzAABFaGdIg4AHJEABW0EY-lKbwxoVoMACiw8pW4pDynEQAHpuZEYHxqOQMLQiNzSOF4JxudYZAKIAAZAAqAFFuTJaJxSNzGf5ucUAEzggCMADYAAwADzN5ow5NQNLwDnpTlm8xZzrAVlkQ0SxjMtJdDk4JnoboWuOsuW4RBKTqDVyIRHIMQ8qEKSvCDkTyeitIcBGo8k4uPIltI4IIFBTpWOmKcV1MtSwe2C9dCKOUaIAkuQYG3205drUJbFTZUke3OEQCLQoNQRrlSI8g6FuNCrrj2ZyeerNqhFJaMMQYIHJ+3aMQQh5LyuwDJ94fa6uHKpB6+J+-OwYe32B4PhxkUcPHHF9p1nedFw5Fdz3rdc4E3VkwG3VAuV5fw5m4UNjyIU88BfBxb2vQir0KDCJVDZ9YNCN9B1o+t6JonoAF9aOYxwBDAZiAF0gA", 
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
# 
# shinyApp(
#   ui = dashboardPage(
#     header = dashboardHeader(userOutput("user")),
#     sidebar = dashboardSidebar(),
#     body = dashboardBody(),
#     title = "User dropdown"
#   ),
#   server = function(input, output) {
#    output$user <- renderUser({
#     dashboardUser(
#        name = "Divad Nojnarg",
#        image = "https://adminlte.io/themes/AdminLTE/dist/img/user2-160x160.jpg",
#        title = "shinydashboardPlus",
#        subtitle = "Author",
#        footer = p("The footer", class = "text-center"),
#        fluidRow(
#         dashboardUserItem(
#          width = 6,
#          socialButton(
#           href = "https://dropbox.com",
#           icon = icon("dropbox")
#          )
#         ),
#         dashboardUserItem(
#          width = 6,
#          socialButton(
#           href = "https://github.com",
#           icon = icon("github")
#          )
#         )
#        )
#       )
#    })
#   }
#  )

