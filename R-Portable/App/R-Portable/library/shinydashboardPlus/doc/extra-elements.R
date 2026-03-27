## ----setup, include=FALSE-----------------------------------------------------
library(bslib)
knitr::opts_chunk$set(echo = TRUE)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0IdmGVoiUpNSUhJOe0ISUkYiamp4oTwHQYBhLp6+gaGCbiJaAjhOdwcCdgBGXAcAJgEk6tCmfP2nKAJiAwaIYJqa7OXB0-P6kjXBgpvQx6ILkgBJchg1w+oVItFI-XuYGsZ2+zxCawc-zgMAcrzwxxunFIUFIZkhcgggUYbwxNWIfTQnDgMmWABUsABVACi72BTkGtN4S1oS04RHgDnIAA9SABCQakhx7QofL4-CBIwGSpyg8FwSHQp6XVGIgHbEky4FYnF4jyDADuLAgfH8BrZTnJ1Ep1OWADFrAAZDIs5XDMCcnkOQN8gXCsUSw01VQfaWymHyoEfO5mlDxuFbO3AuVwxWJ4GqiEpzWw7VbXXI1GZtnG3FLFN8ABmRCrwMdzppHnpzNZbI5XKDvP56rD4rAktjWbTl1zvoL6qLU5I23LKIz6MjmOxtchnBMZ0WnBbHzbqCpHYc7q9Po3NT7geDQ8FcBFo-HGOjoQnKrBhb9xflEZOBOVKMAAbvEywNiYEAEKCJDsHwqDNJsRDNEhpCbFSnCcJcAgOCA+xEPQIHgUy4FkAhEDoQAJMQZATP0jCbARG4mKgcjkP+cIIeeDyLhAaKYXA-SwS6HiIc0tHjL0jEfk4AC+cnDIaREkUEqCMHwpCUTR2aXGsuz7IRxHxKR5HaRJpDUXpLzMRiPBEOaAByRCgg28w4pceY1JwGlaQ2HQAJqoQ4VIieQNI2SETkAA0OAApIeQyWdZ-EGZskqkKwqDzn68DYQEcBHqkJjMHBIQeK8G5KYpKTySIqhgPJAC6QA", 
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
#       radioButtons("controller", "Controller", choices = c(1, 2)),
#       br(),
#       accordion(
#        id = "accordion1",
#         accordionItem(
#           title = "Accordion 1 Item 1",
#           status = "danger",
#           collapsed = TRUE,
#           "This is some text!"
#         ),
#         accordionItem(
#           title = "Accordion 1 Item 2",
#           status = "warning",
#           collapsed = FALSE,
#           "This is some text!"
#         )
#       ),
#       accordion(
#        id = "accordion2",
#         accordionItem(
#           title = "Accordion 2 Item 1",
#           status = "info",
#           collapsed = TRUE,
#           "This is some text!"
#         ),
#         accordionItem(
#           title = "Accordion 2 Item 2",
#           status = "success",
#           collapsed = FALSE,
#           "This is some text!"
#         )
#       )
#     ),
#     title = "Accordion"
#   ),
#   server = function(input, output, session) {
#     observeEvent(input$controller, {
#       updateAccordion(id = "accordion1", selected = input$controller)
#     })
# 
#     observe(print(input$accordion1))
# 
#     observeEvent(input$accordion1, {
#       showNotification(sprintf("You selected accordion N° %s", input$accordion1), type = "message")
#     })
#   }
#  )

## ----message-code, eval=FALSE-------------------------------------------------
# list(
#   author = "David",
#   date = "Now",
#   image = "https://i.pinimg.com/originals/f1/15/df/f115dfc9cab063597b1221d015996b39.jpg",
#   type = "received",
#   text = "A new message"

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0IcAM1NsrCIAd2CamqgCUloSUpNSUhJ7MEY4GCIANzghPAcZrDHJuAd4Tk4A6bAk6vbO7t7+wYhh2JkZ3DmwaxkZFbg1jZntwvaHPZ6IPoGhmZNUOXI50uAFV-lByHcHoEnjsnM9XhATPBGLQCABJCCofrDPhxAAeQJmAFl7utAg5cXA8YggQBGC4wPjuBz0lZQPHMgDM8PaTHysIcfJ5NTM8RJUPubVeDma2VI3GZADYCtKHJxSOCzMyZnIIIFGOcBU5stqwKsyZsVdLRYxxRapaqoP1uERGKbrNQqVAIHE3YZaPECJajUUIR4ZgAmAAMDgAUt6HBHEFGY6gYIaXtLaDANqbuANUJxEAB6YuxRmicgYHrF+Vje7F64VgAyABUAKLFmS0dXF7P+Ys22kAWlpEYAHHix+OMAArVD+DOqpykVioZbhkZwIO0KZnPAhpwzNGcBzyntnsaoajg5ajKDUaisOquuqjOAAfgcre44IA5KeJgQPQcB0HAExQPQnoAIQzCGwrtDadobA60pOvKr6bhkLBQAqfSPkQBAANZLsuAIbpcEacnGCZJlGACsDhpqRqrZrmm75qQhYlmWMgVtQVY1nWqyNnxfBtp23a9v2g6cPEnKjhOU4TnOC4sdKq7rqaclkOprwzAAmkQJiCnAAzxKZYFTBSpCwWA8ECqoNQId0pCeqaNoOMh0L2YUPJyYwUxuh4tRAV0HzsHwWKkBcxlcf0FxyWsHwCA4IA7EQ9ABVM7ZTGQkWYv0AAkozjFMFzpZmTh-ORwIBd5kozOajyzO8JCmqVSxApSHIeFFxU9U5oQAL5DU4mXZXAuUUKQBXRUVpwVQKNW3nVYqkhsnCoTUTUbT5Fwhm1ISbqchJVaExBkDNzJ0Oq227M6mGXAAIlAEzZHp7TkaaAByLSfTUbHkhxBZFqWtAYKgfD9hgxAwMWrq0P4fAPpwxa1LSxa0vRXa1OjtLYzItQEAAnAQkFRoqnL0STADs9BjhGtIyFG2MkyTir0JyJOqYuB7nTUmkUTMow7nuAOhOQeKkO6DgQHAzSQhacECw4Y2hOro07BlWXxDleWzf1pBFSt5BLedptwGttp7Y1qu7RKZ2OuF7WbpbTtZj6VLMkbRU9VaryXeQZA3T2s2HlLMseBq-jNmH90dOgXxHAngOFaQaK3JuozUEo+4B6qN4gdQpoAMJ0MRdx2QdquhKirsUpdwycKwEAEE8NfLmExQRDIpSxIE7CssQufBZcrrej56uvNPzkOedmvqzrk3Tflvs53n5vtDwLR-d0tSouCEWGTE3DV54JjMPsx0shcQumnEtROgJMKZlrhTDY4AhgMNAC6QA", 
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
#     dashboardHeader(),
#     dashboardSidebar(),
#     dashboardBody(
#       fluidRow(
#         actionButton("remove", "Remove message"),
#         actionButton("add", "Add message"),
#         actionButton("update", "Update message")
#       ),
#       numericInput("index", "Message index:", 1, min = 1, max = 3),
#       br(),
#       br(),
#       userMessages(
#         width = 6,
#         status = "danger",
#         id = "message",
#         userMessage(
#           author = "Alexander Pierce",
#           date = "20 Jan 2:00 pm",
#           image = "https://adminlte.io/themes/AdminLTE/dist/img/user1-128x128.jpg",
#           type = "received",
#           "Is this template really for free? That's unbelievable!"
#         ),
#         userMessage(
#           author = "Sarah Bullock",
#           date = "23 Jan 2:05 pm",
#           image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#           type = "sent",
#           "You better believe it!"
#         )
#       )
#     ),
#     title = "user Message"
#   ),
#   server = function(input, output, session) {
#     observeEvent(input$remove, {
#       updateUserMessages("message", action = "remove", index = input$index)
#     })
#     observeEvent(input$add, {
#       updateUserMessages(
#         "message",
#         action = "add",
#         content = list(
#           author = "David",
#           date = "Now",
#           image = "https://i.pinimg.com/originals/f1/15/df/f115dfc9cab063597b1221d015996b39.jpg",
#           type = "received",
#           text = "A new message"
#         )
#       )
#     })
# 
#     output$distPlot <- renderPlot({
#      hist(rnorm(input$obs))
#     })
# 
#     observeEvent(input$update, {
#       updateUserMessages(
#         "message",
#         action = "update",
#         index = input$index,
#         content = list(
#          text = tagList(
#           appButton(
#            inputId = "reload",
#            label = "Click me!",
#            icon = icon("sync"),
#            dashboardBadge(1, color = "orange")
#           )
#          )
#         )
#       )
#     })
# 
#     observeEvent(input$reload, {
#      showNotification("Yeah!", duration = 1, type = "default")
#     })
#   }
#  )

