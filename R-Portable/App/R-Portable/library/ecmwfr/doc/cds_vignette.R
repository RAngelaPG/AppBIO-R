## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# load the library
library(ncdf4)
library(terra)
library(maps)

# grab demo nc files
ncfile <- list.files(system.file(package="ecmwfr"),"*.nc", recursive = TRUE, full.names = TRUE)
ncfile <- ncfile[grepl("cds.nc", ncfile)]

## ----demo request, echo = TRUE------------------------------------------------
# Specify the data set
request <- list(
  "dataset_short_name" = "reanalysis-era5-pressure-levels",
  "product_type"   = "reanalysis",
  "variable"       = "temperature",
  "pressure_level" = "850",
  "year"           = "2000",
  "month"          = "04",
  "day"            = "04",
  "time"           = "00:00",
  "area"           = "70/-20/30/60",
  "format"         = "netcdf",
  "target"         = "era5-demo.nc"
  )

## ----spatial-request, echo = TRUE, eval = FALSE-------------------------------
#  # Start downloading the data, the path of the file
#  # will be returned as a variable (ncfile)
#  ncfile <- wf_request(
#    user = "2088",
#    request = request,
#    transfer = TRUE,
#    path = "~",
#    verbose = FALSE
#    )

## ----spatial-plot, echo = TRUE, figure = TRUE, fig.width = 8, fig.height = 6----
# Open NetCDF file and plot the data
# (trap read error on mac - if gdal netcdf support is missing)
r <- try(terra::rast(ncfile))

if(!inherits(r, "try-error")) {
  terra::plot(r, main = "ERA-5 Reanalysis Demo (2m Temperature 850 hPa)")
  maps::map("world", add = TRUE)
}

