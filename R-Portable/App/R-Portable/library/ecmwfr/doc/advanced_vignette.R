## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

opts <- options(keyring_warn_for_env_fallback = FALSE)

# load the library
library(ncdf4)
library(terra)
library(maps)
library(ecmwfr)

## ----eval = FALSE-------------------------------------------------------------
#  # set a key to the keychain interactively
#  user <- wf_set_key(service = "webapi")

## ----eval = FALSE-------------------------------------------------------------
#  list(stream  = "oper",
#       levtype = "sfc",
#       param   = "167.128",
#       dataset = "interim",
#       step    = "0",
#       grid    = "0.75/0.75",
#       time    = "00",
#       date    = "2014-07-01/to/2014-07-02",
#       type    = "an",
#       class   = "ei",
#       area    = "73.5/-27/33/45",
#       format  = "netcdf",
#       target  = "tmp.nc") %>%
#    wf_request(user = user, path = "~")

## ----eval = FALSE-------------------------------------------------------------
#  # this is an example of a request
#  dynamic_request <- wf_archetype(
#    request = list(
#    "dataset_short_name" = "reanalysis-era5-pressure-levels",
#    "product_type"   = "reanalysis",
#    "variable"       = "temperature",
#    "pressure_level" = "850",
#    "year"           = "2000",
#    "month"          = "04",
#    "day"            = "04",
#    "time"           = "00:00",
#    "area"           = "70/-20/30/60",
#    "format"         = "netcdf",
#    "target"         = "era5-demo.nc"
#    ),
#    dynamic_fields = c("area","day"))
#  
#  # change the day of the month
#  dynamic_request(day = "01")

## ----eval = FALSE-------------------------------------------------------------
#  # creating a list of requests using wf_archetype()
#  # setting the day value
#  batch_request <- list(
#    dynamic_request(day = "01"),
#    dynamic_request(day = "02")
#  )
#  
#  # submit a batch job using 2 workers
#  # one for each in the list (the number of workers
#  # can't exceed 20)
#  wf_request_batch(
#    batch_request,
#    workers = 2,
#    user = user
#    )

