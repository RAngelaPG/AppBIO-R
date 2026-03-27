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
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSdoACADwBaK7u6LlhqAHM47EVaulapNRwVgC8VkJgAEJEAB5WAIKGAJIRuD6O8s76MgAScFAycIzsAmkQvk4uBgDKtIX0LCVlFZlVMtEyHOm+pB6cACScpKxB3mCKnVYgVg0EANbujEQmEDJ2xNREjIhW7txEQwDuvOQAvhGl3VYAZqZ1WESH3uW+r1AE-iSRJqSkJGN-dzuIIAfUUMVS4TAABUiECglZohCwJcXq8rO9PhBvr9-hFGHAYEQAG5wMGxSERLCEknBJGQgjUeScUJQ+juOxyCCeRgXZrojEfWhfH5-CBjAlDLZk8GUsDUqUExEUvBWRnM1kRdl2TgmAgEOCcTh8q6+THC7GivFgEyoOTkcnI3BQgCqdqg5GVTrVTKNmvGHNQjFoMAkFyuqPRTCaV3BzwFfgCCLC5BipAA8j9UD8xuCQUMPXATWjXhEkTMiJ1UqarHV-TBWLK8DWNkzUJwGMmrNCsC6AKL89GMg5Qehdnv9wevVCbDNZnMRGdEUjhkuqXyqVQiThFUmMWwOa4rIX-PjZ0jO5akc-OndGi0CKbpK-n-pL0gHqwE1ZFYzL9ggFcBIAI7sAAhGePz9A24L9K2aA7ioVzvuwjAQFsMDsAATAADDhAjrlYpxbi8L5QXmBaevYX4UIUjDQnAaYAVclHBNRtDXFY7CQaQ0GNrEcFENQbaIY+ETwe2cAyBEVhwNQO5QseQkiVJETIfI5BjAAsqwFZxLQxqqqxhHESIz70DujCkn2pJkNxECvoCwIyrEzqASWtr2nASJjDBKrOuaJD+k5QSrr4pkQOZlnWbZpD2a+BJEqSjpuVcnmFj5ER+d6gXlGE+I0qSYVESRvhEBZu5wDZFBxTx-SSn8BIpU+HnuuQmVgNlkK5f6DXSsVEVRZV1V2XV6UOuCqWtV5mUlr4WX8U6NY9flNptUWzZzVYRCoFiLJhHQQzxgmiaBMEYTcFhYwAHJwIcp2haqlTZAAMqOcnsAAjLevSkGY-pBiGYYopGJ2Uf9q1cjy1Zbb4nBCXUeQFEUrITgONa+IcdSkNwrIACw1oRG7pINpEVVZVWxfFUHZf0xIGZ2cBTUOykIVJn4cVxdW0xJolQrzqlgLJ8nBBESnCWz0lgFc9MdmObEOJz1O8bTsuM2JYBq-LMlyQpES8DIhSRdLJbwEaHgK1YqAaV4EQ6XptaGc6AsyM6ERQKskJa0ERNWDwjzXcuHG0AQHoWuwZucBbzrDKg51QocLAQHw7iQjIJjMFirKfSZqjnJFEBSPw8ToOwli3pVjACGApwALpAA", 
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
# ui <- dashboardPage(
#   title = "Box API",
#   dashboardHeader(),
#   dashboardSidebar(),
#   dashboardBody(
#     tags$style("body { background-color: ghostwhite}"),
#     fluidRow(
#       actionButton("toggle_box", "Toggle Box"),
#       actionButton("remove_box", "Remove Box", class = "bg-danger"),
#       actionButton("restore_box", "Restore Box", class = "bg-success"),
#       actionButton("update_box", "Update Box", class = "bg-primary")
#     ),
#     br(),
#     box(
#       title = textOutput("box_state"),
#       "Box body",
#       id = "mybox",
#       collapsible = TRUE,
#       closable = TRUE,
#       plotOutput("plot")
#     )
#   )
# )
# 
# server <- function(input, output, session) {
#   output$plot <- renderPlot({
#     req(!input$mybox$collapsed)
#     plot(rnorm(200))
#   })
# 
#   output$box_state <- renderText({
#     state <- if (input$mybox$collapsed) "collapsed" else "uncollapsed"
#     paste("My box is", state)
#   })
# 
#   observeEvent(input$toggle_box, {
#     updateBox("mybox", action = "toggle")
#   })
# 
#   observeEvent(input$remove_box, {
#     updateBox("mybox", action = "remove")
#   })
# 
#   observeEvent(input$restore_box, {
#     updateBox("mybox", action = "restore")
#   })
# 
#   observeEvent(input$update_box, {
#     updateBox(
#       "mybox",
#       action = "update",
#       options = list(
#         title = h2("New title", dashboardLabel(1, status = "primary")),
#         status = "danger",
#         solidHeader = TRUE,
#         width = 4
#       )
#     )
#   })
# 
#   observeEvent(input$mybox$visible, {
#     collapsed <- if (input$mybox$collapsed) "collapsed" else "uncollapsed"
#     visible <- if (input$mybox$visible) "visible" else "hidden"
#     message <- paste("My box is", collapsed, "and", visible)
#     showNotification(message, type = "warning", duration = 1)
#   })
# }
# 
# shinyApp(ui, server)

## ----boxTools, echo=FALSE, fig.cap='Box Tools. DFrom left to right: boxLabel, boxDropdown, collapsible and closable buttons, boxSidebar trigger.', fig.align = 'center', out.width='50%'----
knitr::include_graphics("figures/boxTools.png")

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJyduOCgZOEZ3T3lvfRkACQioxnYBXEdQxRlWWK8fAwAhIjzg0JyiAA8KytDSWlJqOFihMABVVDlyB0Vqh05aKPoWdqyQ+qcCaiJOKHoW2IAVLA6AUQmppwB3YdJuWIBGACYt7fDaf25SNrAAVgAGR9Rq8eypziI6ZNTo2IAYtYADIAZU2H3qxGo1DQQ0WrQ8qw25ymUAIjRIxRMpFIJHsYBM3Sg5HGDnayyI-n8SwIykGwzgo0Y7UykMqQxGLFi-VBjOZdW2DmGdxgrDpBk5TLGeHZ9T2MgOsRO91R204P2iAEkIKgcYKheSwER6Jx3pNDUaAHImGD0f5EABmDhNnGiADcSbQSJxEObLU4YHxYo81UKYFABh4js9Q3Kpp7TIiHE9HvHQqohWyLfVULNSAB5HF60gEmS0TikYxEUisuWZjNqqXMgrxIoyPlc9IN7NON2Md3-DyOkwQDHeiDsPgl3Au4s42duzhDEgCBwgSGuj1BVCMPil6c4gAkYolMmbLAEFhzRHnpCP5cr1duAB4ALQORgUNLP9gbnNhBWpaMBARCMDAU66serpXpCAC+16VFuA5wOsg5kJBJZHkSPRwLO-5TDhJJwKU1SdtK6TtKeygXiyYANk4CHZHBIgCGAcEALpAA", 
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
#     body = dashboardBody(
#       box(
#         title = "Update box sidebar",
#         closable = TRUE,
#         width = 12,
#         height = "500px",
#         solidHeader = FALSE,
#         collapsible = TRUE,
#         actionButton("update", "Toggle card sidebar"),
#         sidebar = boxSidebar(
#           id = "mycardsidebar",
#           width = 25,
#           sliderInput(
#             "obs",
#             "Number of observations:",
#             min = 0,
#             max = 1000,
#             value = 500
#           )
#         ),
#         plotOutput("distPlot")
#       )
#     ),
#     sidebar = dashboardSidebar()
#   ),
#   server = function(input, output, session) {
#     observe(print(input$mycardsidebar))
# 
#     output$distPlot <- renderPlot({
#       hist(rnorm(input$obs))
#     })
# 
#     observeEvent(input$update, {
#       updateBoxSidebar("mycardsidebar")
#     })
#   }
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0IdFAA9gmprSWlJqOHcHITAAYWoiTih6dodyuocAd1buT0YiVBkiCYhugqbQgn7B4Y6PABUsAFUAUTX1yezSGY8ARgAmM-XOUihSM07uiZYIPn9V6qanCIdBk0Vi8U6ADFrAAZDKnAE1YjUahoTgMEb7I4IwrrGRzBZLCAAWQoJk69QAIgTFstGucnFSaUSAJLkGD2Xp0AgAawc8FWDmyHzA+PmtIgbLgMEFtGIIQ8cpInO4MUYpG6AiSiKaTPFrPZnNa0ocd0F3EYcAAZiLuKRSKhOIgAPTOibujD+IhEfztDDEGDOzWPc5iwnLSm0ABu2Xi+R1NT14clhu6xpgDgAzLL5Z0lRBOVdNaoGdrcU1usTWLUiHVujqS6FG2WnJx4lGIR4rSYIAQWsq+KgTKRcA4iMOh6QBA4QACiPQ24wO8cO2R2IPhwASMMSqUwUez8ut7hLAByRBaVrlr1oyu60RRREFMhMzH7CocN1HpFYqF2XTAeBOEGQJNQBABfRtwJEAQwHAgBdIA", 
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
#       box(
#         title = "Closable Box with dropdown",
#         closable = TRUE,
#         width = 12,
#         status = "warning",
#         solidHeader = FALSE,
#         collapsible = TRUE,
#         dropdownMenu = boxDropdown(
#           boxDropdownItem("Click me", id = "dropdownItem", icon = icon("heart")),
#           boxDropdownItem("item 2", href = "https://www.google.com/"),
#           dropdownDivider(),
#           boxDropdownItem("item 3", icon = icon("th"))
#         ),
#         "My box"
#       )
#     )
#   ),
#   server = function(input, output) {
#     observeEvent(input$dropdownItem, {
#       showNotification("Hello", duration = 1, type = "message")
#     })
#   }
# )

## ----user-description, eval=FALSE---------------------------------------------
# userDescription(
#   title = "Nadia Carmichael",
#   subtitle = "lead Developer",
#   type = 2,
#   image = "https://adminlte.io/themes/AdminLTE/dist/img/user7-128x128.jpg",
# )

## ----userBox, echo=FALSE, fig.cap='Some userBox components', fig.align = 'center', out.width='100%'----
knitr::include_graphics("figures/userBox.png")

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0OdOePKAD2CamtJaUmo4d1r4gBE4TgJGWlQ2kmaW1vbO7qEwADlY2igHAGEWGFoCbig4ajmCyZrOE3o2jq6POc7Yh36ANz2iVHiD6qPSVhfugCZDo6ctBgAUuDjm3FIpFQnEQAHpYbFNqJyBhaERYaRuHB4JxYdYZEiADIAFQAorCZLROKRYUD-LCzPEAOwAWgAjD8ABwNDmcjAAK1Q-jehUmSXex1IUFIZlmYAA7iwIHxhXgJaFoPdDLRqNROBMAdk5agddQAJLkGAigFOTXa3WcC3Yg02zoAM1IcqdMAcbOtNqcxGoREYcv8jDgFH9AeG-gh3TZAAZ1TVxaKjnbTY7LS6Ae7PVcwN6HD9ozagyG5RGZGXDcQQh4tuM5lAIP5OiyZER5RA5mmA05Y-HC0mAKRzFOhVQA-tHN1EIjkUOF4lYhzzxfxBxYiMAQgn6ac06n-yc5xmhcZjEaB6P-zqjEey-XJggBDGEHYfFQJlIuAcRC-j+pACA4IDvEQ9APo8pKPGQX4QMBAAkJq6t6-7gYeDi0G6DgIchqHmpa7geD8oGYQCPDdgAshUUDULmLQwHR1C9Mswb+PYYDWA4zFyPsYDHmK6oAL7vCJFjppB0FwLBFCkPhv4oaa6FgeqVHyvMi44Vs0pop+k6cKgwxkG6XEAJqAQ4BB0AQADWcAyABISEQ48wAA0OKOnAHNhiFKYR3qzpMnzfIWiqMMqba1qkJjMB+CbqkJDgSSkYkQAIYAiQAukAA", 
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
#       userBox(
#         title = userDescription(
#           title = "Nadia Carmichael",
#           subtitle = "lead Developer",
#           type = 2,
#           image = "https://adminlte.io/themes/AdminLTE/dist/img/user7-128x128.jpg",
#         ),
#         status = "warning",
#         navPills(
#           id = "pillItem",
#           navPillsItem(
#             left = "Item 1",
#             color = "green",
#             right = 10
#           ),
#           navPillsItem(
#             left = "Item 2",
#             color = "red",
#             icon = icon("angle-down"),
#             right = "10%"
#           )
#         ),
#         footer = "The footer here!"
#       )
#     ),
#     title = "userBox"
#   ),
#   server = function(input, output) {
#     observeEvent(input$pillItem, {
#       if (input$pillItem == 2) {
#         showModal(
#           modalDialog("A modal")
#         )
#       }
#     })
# 
#     observeEvent(input$pillItem, {
#       showNotification(
#         sprintf("You clicked on pill N° %s", input$pillItem),
#         type = "warning",
#         duration = 1
#       )
#     })
#   }
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0IcoAlJaElKTUlISezATVDlyAH1FAA8hPAdhgFVuqHJhpOrQs3jygeCamuz3Uc7OeMHhgtXQhtJqOA2FxgAROE4CRlpUBva5g6OTjeGAOVjaKAcAYRYMFoBG4UDg1D2z1WnBM9Fepw8wxOsQcVwAbuCiKh4pDCgcnKRWNiNgAmfb4py0GABBGbbitVCcRAAemZsSBonIGEazNI3Dg8E4zOsMg5ABkACoAUWZMlonFIzKp-mZ5wA7ABaACMJIAHAMdbqMAArVD+XEU2Z46GkKZmd5gVB3alsC34-zMOUUUgbCVYMZS8kHXIEADWHqIJggMgddH89I19FMcDdwaIAwlRCI1CyAC9acMBhC8FCnMMMkR4A5yAMffzGHAAISp1YAMyz5EYDol-Ic7aInYc9abwyhqhqVpq8Id5yWo8Kk4c20YGK7HlbUfqjQg7D4qBauAckdI+9IAgcIDmRHoy4xUoxZF3EFPABIuj04P104fL9anO+pjgJYVnxYZzl2Es-1COpHhCRFOkmaZIIpLFYM4DY6AVECKWnDxziuG47gebdsIpataGOAswAAKRiEI52QsinBhOEKLeeD0UxbFGBbHCiVpLUgzIqkaQdekTyZVl2T4aguR5PkBWuYVRT4SUZTlBUlRgFVzgAZm1PUDT1E0zV41ZFwpBU7XQ+CGxkMyahDcNGEjaMNg+MYxTFISDgAd2yPkNgAFlLBxx3MuYAF9wsikQBDASKAF0gA", 
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
#       actionButton("update_box", "Update"),
#       userBox(
#         id = "userbox",
#         title = userDescription(
#           title = "Nadia Carmichael",
#           subtitle = "lead Developer",
#           type = 2,
#           image = "https://adminlte.io/themes/AdminLTE/dist/img/user7-128x128.jpg",
#         ),
#         status = "primary",
#         gradient = TRUE,
#         background = "light-blue",
#         boxToolSize = "xl",
#         "Some text here!",
#         footer = "The footer here!"
#       )
#     ),
#     title = "userBox"
#   ),
#   server = function(input, output) {
#     observeEvent(input$update_box, {
#       updateBox(
#         "userbox",
#         action = "update",
#         options = list(
#           title = userDescription(
#             title = "Jean Box",
#             subtitle = "Developer",
#             type = 1,
#             image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#           ),
#           status = "red",
#           background = NULL,
#           width = 4
#         )
#       )
#     })
#   }
# )

## ----user-block, eval=FALSE---------------------------------------------------
# userBlock(
#   image = "https://adminlte.io/themes/AdminLTE/dist/img/user4-128x128.jpg",
#   title = "Social Box",
#   subtitle = "example-01.05.2018"
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0IdOIgJaKGpygA9gmprs9wchMDqGpsUW3oKO0NJaUmo4brN40up6gGt2sdDaGACZj17uUlJUTkQAemPYmD5qcgxaImPSbjh4TmPrGQuIABkAFQBRY5ktE4pGOG38xzmjAALABaACMACYABwtRFIjAAK1Q-hG1TWEym2x6YAy9Ua1AcrVxhTWtRM9AJ026vTgLVgqGmMIADHCMFyAKwYBE8pG9PE1JLi0JQAgTEilEz7Ej2MAmVBycgAfSGI2JWDgADNGHAeL1JTSOr1SfAHOQWqQHI9jQBCam0pj5KVOUhQeicOCkPwQODUVa0n30IMhsO04kAOTgAHdOG7Y04oPsZdx4GQFssY2mNltmWA9gcjqdE1WMJxYDAiPRaNMMMQYMdE6gYcQyBQQWrFrEXsK4QB2Y5wuHHADCsBDGFQEBxeC9+MmTJ2YG+JtIqbTTm4xoNJbLhxOx38RCI-mbrd3e96314nAcQNtjwc3fIZDFFtpqjT5pphGUahiuNS9AAsiataBCmy6-mskJQZwMEmgWsaJtkDzdIiox7k4wIZmYJZyIu8R3mmSHQVs6FplAircEQjAltY0xshAcTMYYtDxAQcAUXuGpEr0woOAAUlAIQIogXJcg4qAwAJhabIEx77KepznJc1y3PcjzPK87x8D8-yAsCoIwOCkJwvCyKosimLYkpsakKwqDCWAxp8bQABucAyM5tK9AAks+DyvuQMAchmMzGk01CsA4BpMUlxpwAA-A4j4ZgA5M+JgQPQIY8T5vrTK6YBgR0gH4c4-qMMhqG0bG9EPClG4ZCwUDcJSJjUIsBBLIFaxCSWCIAMziZJDjSQK8mKfBtXrCpHknhWZxGaIOl3A8TwmoZHwmQCQIgmCEL1eNtkomijlLnhS2ue5Jb+t+i1LcSACaRAmA4RX7PEv3FXAfkvqQFVVRKEMOP+sYw2MNVjNQaAchwAAaOGIHCXK4A4ABiACqcbdAaBWyrcEDsLQAgOCAVVDFORAwDmpDNZ0K1qeWZ5aVtcA3Dt+n7W8h1-Md5lnZCl1ovZ6JYndUOMkSqDyOQKoM0zva6lT92xqNG48nygrDqKb1pkrwJwFyKqPjMms470MLYa2zNmlVcMdAAvgjHTJUQ5DMRu1tJZefuOvEcDgwhbtewrJakgMFJUpVhQ1fVfn+0lpNyhTfCoIqOPfQcirU7TFoNqncC-H5ZCUxAuekAAJGqQnakQLQ4yXYxNzFrTNb0-TkjqJsdDKWcll35DOUQqBZ8+Hh0MCrO2muRKQnmg2L04RaqRua1c5tVy87pu0GULxki2Zp2Wed8SS3ZN2y8N3rL7HZJNJSrfOOqMUBUPtKcPSMcNykAqFAVgP4AIrjdk4N2nsUjuxEAIMA7sAC6QA", 
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
#       socialBox(
#         id = "socialbox",
#         title = userBlock(
#           image = "https://adminlte.io/themes/AdminLTE/dist/img/user4-128x128.jpg",
#           title = "Social Box",
#           subtitle = "example-01.05.2018"
#         ),
#         actionButton("update_box", "Refresh"),
#         "Some text here!",
#         br(),
#         tabsetPanel(
#           tabPanel(
#             "News",
#             attachmentBlock(
#               image = "https://www.sammobile.com/wp-content/uploads/2017/11/Camel.png",
#               title = "Test",
#               href = "https://google.com",
#               "This is the content"
#             )
#           ),
#           tabPanel(
#             "Messages",
#             userMessages(
#               width = 12,
#               status = "danger",
#               userMessage(
#                 author = "Alexander Pierce",
#                 date = "20 Jan 2:00 pm",
#                 image = "https://adminlte.io/themes/AdminLTE/dist/img/user1-128x128.jpg",
#                 type = "received",
#                 "Is this template really for free? That's unbelievable!"
#               ),
#               userMessage(
#                 author = "Sarah Bullock",
#                 date = "23 Jan 2:05 pm",
#                 image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#                 type = "sent",
#                 "You better believe it!"
#               )
#             )
#           )
#         ),
#         lapply(X = 1:10, FUN = function(i) {
#           boxComment(
#             image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#             title = paste("Comment", i),
#             date = "01.05.2018",
#             paste0("The ", i, "-th comment")
#           )
#         }),
#         footer = "The footer here!"
#       )
#     ),
#     title = "Social Box"
#   ),
#   server = function(input, output) {
#     observeEvent(input$update_box, {
#       updateBox(
#         "socialbox",
#         action = "update",
#         options = list(
#           title = userBlock(
#             image = "https://adminlte.io/themes/AdminLTE/dist/img/user3-128x128.jpg",
#             title = "Social Box updated",
#             subtitle = "today"
#           )
#         )
#       )
#     })
#   }
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABA5O0HAXge7ui5YagBzOHsIJycvHwMACTgoGThGdgFcR1Dw-RkAZVo4+hZE5JDU+W90gCEiGQ4U0IcAM1NsrCIAd2CamuJTGAg29prm7NJudwcANgK+mpcAeRNSVDn7FAJSWgA3OAB9Tmy4ITAk6smoFdoSUrnSEiWr-39qPbwHfYAVIjuHurpUB0UAD32h0Kky+tFQ5T+vRBTmyI32MFY9TB-32E2hTlIjFod3icLA3CIG0YqKOIIGMiGIwAjAAmNHo2qMEikEYydZQ9EOAjUeScPHkP6kAC0BAo5GJeFJ0O4VKWADFvg4SA4CUTAfTObQYP4uIwCHjuKR5pxEAB6U2oHmignUGQYYim0YABidfwALC7AVK+kDOb8TgBrVns72Tbm8-lwQUisXxEnAzncODYw14gDMLtQAMlCfR5MpHn2VJdAFJ436nDL5Yrlaq4wcNejUEsALJERhwBx8WrtmBQVYkDBDjBe3OTVToifjxuhTombocvr54YecahpwzOYLUhLE6rDbbXabGnq9eg8FESFnmEyPEIpGof4nnN+5fUunXupMsjBtaL6HhpwfKFmAArCqKZD1jOIJVvsCpgkqITcrQBABqeY4glqOqcHqBpGqgJrmpaJxwDadoOs6roek6o5+r6nK5Khv7-iCgHAc8oFRuBsYStBkxJimLIgRmTpZuWFZviBxZOmWL4VrBYDwT8yrIah6EVg4zb7G2HZdhAPaMH2A4QMOI4HJ+U4gpZ7TWU41lTvRnDxESIy1CYECnNcfDbrgSpbnMvlOUBZwQAIDggEcRD+aQAAke7rFsOxxA4AA8QoOB2EBxIwACqACS7ARWOSWdmlXa1A47DeXMMX3t8-xhfsjLMvsDhwNQTkcYxaFgFKaTKKUsSBFwuy+Z07Z4vQpiPNZAC+FgJlF8w1fFB4lceqXpZl2X5YVUolZt5WVdVsV1cil40o1YDNWQrXtZ1+zdfsfXFBEMiDTIw0lWNRDUBNIFTSYM1HPNIiRfQTmMBsACiGxkFVEDbjFtz3HAvlFX0JioHI5BKRCSxnY+l60aEoOFLNIgCGAs0ALpAA", 
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
#         column(
#           width = 6,
#           uiOutput("active_side"),
#           actionButton("toggle", "Toggle flip box"),
#           flipBox(
#             id = "myflipbox",
#             trigger = "hover",
#             width = 12,
#             front = div(
#               class = "text-center",
#               h1("Flip on hover"),
#               img(src = "https://placehold.co/600x400")
#             ),
#             back = div(
#               class = "text-center",
#               height = "300px",
#               width = "100%",
#               h1("Flip on hover"),
#               p("More information....")
#             )
#           )
#         ),
#         column(
#           width = 6,
#           uiOutput("active_side_2"),
#           flipBox(
#             id = "myflipbox2",
#             width = 12,
#             front = div(
#               class = "text-center",
#               h1("Flip on click"),
#               img(src = "https://placehold.co/600x400")
#             ),
#             back = div(
#               class = "text-center",
#               height = "300px",
#               width = "100%",
#               h1("Flip on click"),
#               p("More information....")
#             )
#           )
#         )
#       )
#     )
#   ),
#   server = function(input, output, session) {
#     output$active_side <- renderUI({
#       side <- if (input$myflipbox) "front" else "back"
#       dashboardBadge(side, color = "blue")
#     })
# 
#     output$active_side_2 <- renderUI({
#       side <- if (input$myflipbox2) "front" else "back"
#       dashboardBadge(side, color = "blue")
#     })
# 
#     observeEvent(input$toggle, {
#       updateFlipBox("myflipbox")
#     })
#   }
# )

