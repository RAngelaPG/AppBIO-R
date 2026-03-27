## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(attachment)

## ----echo=FALSE---------------------------------------------------------------
attachment:::extra_dev_pkg

## ----eval=FALSE---------------------------------------------------------------
# create_renv_for_dev() # with all default above
# create_renv_for_dev(dev_pkg = "attachment") # with {attachment} only
# create_renv_for_dev(dev_pkg = c("_default", "DT")) # for all default and {DT}

## ----eval=FALSE---------------------------------------------------------------
# create_renv_for_prod(output = "renv.lock.prod")

