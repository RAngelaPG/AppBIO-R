## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

# load the library
library(ncdf4)

# disable render of chunks on mac
if(tolower(Sys.info()[["sysname"]]) == "mac") {
   knitr::opts_chunk$set(eval=FALSE)
}

## ----echo = FALSE-------------------------------------------------------------

# grab demo nc files
ncfile <- list.files(system.file(package="ecmwfr"),"*.nc", recursive = TRUE, full.names = TRUE)
ncfile <- ncfile[grepl("cds_workflow.nc", ncfile)]

# open the netcdf file and print the meta-data
f <- ncdf4::nc_open(ncfile)

# read in the temperature data stored in field "t2m"
# and the dates from the "time" field
temp <- ncdf4::ncvar_get(nc = f, "t2m")
time <- ncdf4::ncvar_get(nc = f, "time")

# get the starting point of the time series
# and add the increments (time)
start <- ncdf4::ncatt_get(f, "time")$units
start <- as.Date(start, format = "days since %Y-%m-%d 00:00:00")
time <- start + time

# close the file once done
ncdf4::nc_close(f)

# plot data
plot(time, temp, ylab = "Temperature (K)")

