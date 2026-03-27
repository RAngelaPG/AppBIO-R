## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(attachment)

## ----examples-create_dependencies_file----------------------------------------
# Create a fake package
tmpdir <- tempfile(pattern = "depsfile")
dir.create(tmpdir)
file.copy(system.file("dummypackage",package = "attachment"), tmpdir,
          recursive = TRUE)
dummypackage <- file.path(tmpdir, "dummypackage")

# Create the dependencies commands but no file
create_dependencies_file(
  path = file.path(dummypackage,"DESCRIPTION"),
  to = NULL,
  open_file = FALSE)

# Create the dependencies files in the package
create_dependencies_file(
  path = file.path(dummypackage,"DESCRIPTION"),
  to = file.path(dummypackage, "inst/dependencies.R"),
  open_file = FALSE)
list.files(file.path(dummypackage, "inst"))
# browseURL(dummypackage)

# Clean temp files after this example
unlink(tmpdir, recursive = TRUE)

