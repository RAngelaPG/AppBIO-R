## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  eval = FALSE
)

## ----setup, eval=TRUE---------------------------------------------------------
library(gfonts)

## -----------------------------------------------------------------------------
#  setup_font(
#    id = "roboto",
#    output_dir = "www",
#    variants = "regular"
#  )

## ---- eval=gfonts:::is_service_ok()-------------------------------------------
head(get_all_fonts())

## -----------------------------------------------------------------------------
#  use_font("roboto", "www/css/roboto.css")

## -----------------------------------------------------------------------------
#  fluidPage(
#  
#    tags$link(rel = "stylesheet", type = "text/css", href = "css/roboto.css"),
#    tags$style("body {font-family: 'Roboto', sans-serif;}")
#  
#  )

## -----------------------------------------------------------------------------
#  download_font(
#    id = "roboto",
#    output_dir = "azerty",
#    variants = c("regular", "300italic", "700")
#  )

## -----------------------------------------------------------------------------
#  generate_css("roboto", "regular", font_dir = "path/to/font")

## ---- eval=gfonts:::is_service_ok(), echo=FALSE-------------------------------
cat(generate_css("roboto", "regular", font_dir = "path/to/"))

