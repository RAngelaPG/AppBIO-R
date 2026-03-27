## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(dockerfiler)

## ----example-dock_from_renv---------------------------------------------------
#' \dontrun{
#' dock <- dock_from_renv("renv.lock", distro = "xenial")
#' dock$write("Dockerfile")
#' }

