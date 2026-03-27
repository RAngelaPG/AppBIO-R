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

# grab demo nc files
ncfile <- list.files(system.file(package="ecmwfr"),"*.nc", recursive = TRUE, full.names = TRUE)
ncfile <- ncfile[grepl("webapi.nc", ncfile)]

## ----eval = FALSE-------------------------------------------------------------
#  # set a key to the keychain
#  wf_set_key(user = "khrdev@outlook.com",
#             key = "XXXXXXXXXXXXXXXXXXXXXX",
#             service = "webapi")
#  
#  # you can retrieve the key using
#  wf_get_key(user = "khrdev@outlook.com")
#  
#  # the output should be the key you provided
#  # in this case represented by the fake X string.
#  # "XXXXXXXXXXXXXXXXXXXXXX"

## ----eval = FALSE-------------------------------------------------------------
#  # set a key to the keychain
#  wf_set_key(service = "webapi")

## ----eval = FALSE, message=FALSE, warning=FALSE-------------------------------
#  # this is an example of a request
#  request <- list(stream  = "oper",
#                     levtype = "sfc",
#                     param   = "167.128",
#                     dataset = "interim",
#                     step    = "0",
#                     grid    = "0.75/0.75",
#                     time    = "00",
#                     date    = "2014-07-01/to/2014-07-02",
#                     type    = "an",
#                     class   = "ei",
#                     area    = "73.5/-27/33/45",
#                     format  = "netcdf",
#                     target  = "tmp.nc")
#  
#  # an example download using fw_request()
#  # using the above request list()
#  ncfile <- wf_request(user    = "khrdev@outlook.com",
#             request  = request,
#             transfer = TRUE,
#             path     = tempdir(),
#             verbose  = FALSE)

## ----fig.width = 7, fig.height = 7--------------------------------------------
# Open NetCDF file and plot the data
# (trap read error on mac - if gdal netcdf support is missing)
r <- try(terra::rast(ncfile))

if(!inherits(r, "try-error")) {
  terra::plot(r[[1]], main = "2 meter temperature (Kelvin)")
  maps::map("world", add = TRUE)
}

## ----mars example, eval = FALSE-----------------------------------------------
#  # this is an example of a request
#  request <- list("dataset" = "mars",
#                     "class"   = "od",
#                     "date"    = "2019-01-22",
#                     "expver"  = "1",
#                     "levtype" = "sfc",
#                     "param"   = "167.128",
#                     "step"    = "0",
#                     "stream"  = "oper",
#                     "time"    = "12",
#                     "type"    = "fc",
#                     "grid"    = "0.125/0.125",
#                     "area"    = "50/0/30/15",
#                     "format"  = "netcdf",
#                     "target"  = "mars-data.nc")
#  
#  # Submit mars request and wait for netcdf file
#  wf_request(user    = "khrdev@outlook.com",
#             transfer = TRUE,
#             path     = "~",
#             request  = request,
#             verbose  = FALSE)

## ----eval = FALSE-------------------------------------------------------------
#  # user info
#  user_info <- wf_user_info(user = "khrdev@outlook.com")
#  head(user_info)
#  
#  # services
#  services <- wf_services(user = "khrdev@outlook.com")
#  head(services)
#  
#  # datasets
#  datasets <- wf_datasets(user = "khrdev@outlook.com")
#  head(datasets)

