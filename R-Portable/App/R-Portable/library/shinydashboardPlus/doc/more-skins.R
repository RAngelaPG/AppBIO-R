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
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVVP4BBdOxEACRydr2AvPd3dFyw1ADmcHYQjiHccFAycIxuHvJe+jIAEuGRjOwCuA4h9py0kfQsMZ7eBgDKeXAFaRlZIYoyrEVxJTIAQkQN6ZnB2cRkjETUVU16ygDCJKQDQyzsxNTUaJxwMjEAYlYAMqUAorg5ANZ8pXDUcASkRNU1PSGktKRnMUJgpUfBJ2cXVy+1NyHLRgANyiMQAZiYIBdaCR2HxUCZSPsiIiEaQBPYQPYAL4OARgbEAXSAA", 
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
#      header = dashboardHeader(),
#      sidebar = dashboardSidebar(),
#      body = dashboardBody(),
#      controlbar = dashboardControlbar(collapsed = FALSE, skinSelector()),
#      title = "Skin Selector"
#    ),
#    server = function(input, output) { }
#  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABE5O0HAXge7ui5YagBzOHsIJ1DOAGs+dwchMBhaGQhaf25SWNxHUIduOCgZOEZorx8DAAlc-MZ2AQyQrM4EuHoWIvlvfRkAZUbmqprM0MUZVla9ZQAhImHgrMGiAA8Z2dDSWlJqOGjYgFkRxXnYgazVWf660OIyRiJqXtH25QBhElJr25Zq2tmAMyIickKHmKHQAYn8AZ8jg5VutNh4dgkkilSA5OpEIIdzmcwgUAG4FaLfEwQAirEjsPioEykXAOIjUqmkAQOEAOAC+jgEYDZAF0gA", 
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
#    ui = dashboardPage(
#      skin = "midnight",
#      header = dashboardHeader(),
#      sidebar = dashboardSidebar(),
#      body = dashboardBody(
#        box(
#          title = "My box"
#        )
#      ),
#      controlbar = dashboardControlbar(),
#      footer = dashboardFooter(),
#      title = "Midnight Skin"
#    ),
#    server = function(input, output) { }
#  )

## ----midnight-skin, echo=FALSE, fig.cap='Midnight skin overview', fig.align = 'center', out.width='100%'----
knitr::include_graphics("figures/skin-midnight-overview.png")

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVSKn8AgunYiABE5O0HAXge7ui5YagBzOHsIJ1CYGXcHABUsAFUAUVxHUIdOAGs+SKEwelM4bKSQlO44KBk4RkivHwMACVLyxnYBQpTU2nL6Fir5b30ZAGUOuC6mluTQxRlWHr1lACEiaeC2p0UAD3ZSWlJqOCywAFkZjezVNvGi0OIyRiJqUdm+5QBhElI7h5Zm1pSAMyIRHIlQ81X6ADFAcCfhMnNtdvsPNlDlBgbQoNRPHBOLR-BBshNLqFOBUAG4VSJ-EwQAjbEjsPioEykXAOIjMpmkAQOEAOAC+jgEYD5AF0gA", 
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
#    ui = dashboardPage(
#      md = TRUE,
#      skin = "blue",
#      header = dashboardHeader(),
#      sidebar = dashboardSidebar(),
#      body = dashboardBody(
#        box(title = "My box")
#      ),
#      controlbar = dashboardControlbar(),
#      footer = dashboardFooter(),
#      title = "Material design"
#    ),
#    server = function(input, output) { }
#  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
card(
  shiny::tags$iframe(
    class = "html-fill-item", 
    src = "https://shinylive.io/r/app/#h=0&code=NobwRAdghgtgpmAXGKAHVA6ASmANGAYwHsIAXOMpMAGwEsAjAJykYE8AKAZwAtaJWAlAB0IdJiw48+rACZQe9IixnDRDZmy69+chUsYyACtQCunVWI0cAZozg9VAYgAEBO1HLPS3OF5-xnAHdab2coV1Z6OEZUEwgAa1ciaiJGZ1QoajhSchFvOACAHgBaV3dyAH18+HYRZzCZGD5qSuIUxlqIevqAczsKZwBeZyEwRwBmaymAJhlR3Dru+lNfYdHHaYA2cYBWKfnF+rsZIZGwIVJHGfpp+gOu7tY4ahTA0-XrOCnNuHvu52sJgI3E4tHCazGUx29GsAA4-t1oAA3VjvMbjADsABYCABOaYI+qxGJZNGOTYARh2BBhhOcMBYRBIZPomx2uLxdLoPW4pAqyxMqzOjh2m1xFIIBFGiwECweUEazUqoJkcHoLE6-zkjHi-J6ZIAIrCDQBRE24una3XcIhI6J6smwikAQQpAGEKZaWLq2qkydMTeMsViAAzSh6yxYKpqiSo9FLq6ia7rEMgUPn0fUQxxQXN58P1VSqEQuNxwDy+UiMBjLXzWP2KAAeznjRETSQg1loPRETYqqa7+pKXhrWUQiCro7gyYAfuqCPE+kQ4jJcM4Z9QoFFqABlUgeMxy+qjY7zM6BFgQPg9vCLUbE1BZM+jThAgj2TgI0Z9OAUZ9gVBqwZNgvzAJ4XiIQJ-zkCAemiUDAWBUEoH-Ph63DYsIBcJsASgAhSFSVFAQgfDaBIXsiEbCprDwgi2GcYdiNIkh2HnRdGGXCBV2cTdtz3A9zGcEBFibZN6lIEJSQhN1ImiWIEmcAAhSjnweeo2k3VBQVrU4ABUsAAVRNNdDmcNilxXU5zI4lcj26Hxu15NFphDENUEbLkt2eKzKIAGS8pMKTXXjnn40gzFUQsRAAXxEEt0gZVB6SgHpaClCAm3sBjSmJRhGHHVBEtYyj+xIQc1z7Gj8MIzCpH4Z10DEkxaFOXRuEUZRDBS6dTIBOweF0-whWqOA7O6Th4j4NEBV+W81PqIhUAkkhOFOOhOFILhaFVdVGBNRsMi4gB5CAAAlbWiPTDJNSN5ucHwFUu4Y2o6gxTvLVUOl6+oZA41AZEgiAAFkKBMMT-nE1hUCFUZ4E4Thus-PARju-51RkOCwrMNFX0lD86QhuGEbggBJcgYHBiH6lsIgYDRHcTHQVJSGcXTyxgAmqaJ7q0UG2hVv5vxfFTcgyGcIhrDCekP26jBOYhiSAghHZ6T4JHvu6W6qe6bnSfJymqZpumIQZpnGBZtnYHl-5dZhsA+YF1b8g7UWWYlsIICIfI0ltuW5u1-5Fbt6Z7uXRh1dRzWNci-4tf+FU1RYVr5Ha-QZB3bbE6+yOE92kGIDBjX6ngAuyYKWowDLunPWR9G4P87dTjSljRge83RllMyFTgt1kj9CFuV5YoZo7saudBquK6r5wCVr7u4Ab7zhmbiAK9QJlR67jG4F79o0WApkIA76PerjpYiBkVFnpT16ZGUy-isbewz-Ukgq2SXbk70ZRe7IDjqF2uwF+I5SBSTOAAMX6twIWMYbyHDPpwaIdo0jDCYstVefBYikDXMuUgWCBBCWcLFLoAgwDRQALpAA", 
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
# library(fresh)
# # create the theme with a cyberpunk color palette
# theme <- create_theme(
#   adminlte_color(
#     green = "#3fff2d",
#     blue = "#2635ff",
#     red = "	#ff2b2b",
#     yellow = "#feff6e",
#     fuchsia = "#ff5bf8",
#     navy = "#374c92",
#     purple = "#615cbf",
#     maroon = "#b659c9",
#     light_blue = "#5691cc"
#   ),
#   adminlte_sidebar(
#     dark_bg = "#D8DEE9",
#     dark_hover_bg = "#81A1C1",
#     dark_color = "#2E3440"
#   ),
#   adminlte_global(
#     content_bg = "#aaaaaa"
#   )
# )
# 
# # create tribble for box global config
# box_config <- tibble::tribble(
#   ~background, ~labelStatus,
#   "red", "warning",
#   "purple", "success",
#   "green", "primary",
#   "yellow", "danger",
#   "fuchsia", "info"
# )
# 
# # box factory function
# box_factory <- function(background, labelStatus) {
#   box(
#     title = "Cyberpunk Box",
#     collapsible = TRUE,
#     background = background,
#     height = "200px",
#     label = boxLabel(1, labelStatus)
#   )
# }
# 
# # pmap magic
# boxes <- purrr::pmap(box_config, box_factory)
# 
# shinyApp(
#    ui = dashboardPage(
#      freshTheme = theme,
#      skin = "blue",
#      options = list(sidebarExpandOnHover = TRUE),
#      header = dashboardHeader(
#        dropdownMenu(
#          type = "messages",
#          badgeStatus = "success",
#          messageItem(
#            from = "Support Team",
#            message = "This is the content of a message.",
#            time = "5 mins"
#          ),
#          messageItem(
#            from = "Support Team",
#            message = "This is the content of another message.",
#            time = "2 hours"
#          )
#        )
#      ),
#      sidebar = dashboardSidebar(
#        sidebarMenu(
#          menuItem("Item 1", badgeLabel = icon("heart"), badgeColor = "light-blue"),
#          menuItem("Item 2", badgeLabel = icon("poo"), badgeColor = "maroon")
#        )
#      ),
#      body = dashboardBody(boxes),
#      controlbar = dashboardControlbar(),
#      title = "Fresh theming"
#    ),
#    server = function(input, output) { }
#  )

