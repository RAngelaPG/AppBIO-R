#' @tags with
#' @tags %dopar%
#' @tags sequential multisession

library(doFuture)

foreach::registerDoSEQ()
stopifnot(foreach::getDoParName() == "doSEQ")

with(registerDoFuture(), {
  stopifnot(foreach::getDoParName() == "doFuture")
})
stopifnot(foreach::getDoParName() == "doSEQ")

with(registerDoFuture("%dopar%"), {
  stopifnot(foreach::getDoParName() == "doFuture")
})
stopifnot(foreach::getDoParName() == "doSEQ")

with(registerDoFuture("%dofuture%"), {
  stopifnot(foreach::getDoParName() == "doFuture2")
})
stopifnot(foreach::getDoParName() == "doSEQ")


local({
  with(registerDoFuture(), local = TRUE)
  stopifnot(foreach::getDoParName() == "doFuture")
})
stopifnot(foreach::getDoParName() == "doSEQ")

local({
  with(registerDoFuture("%dopar%"), local = TRUE)
  stopifnot(foreach::getDoParName() == "doFuture")
})
stopifnot(foreach::getDoParName() == "doSEQ")

local({
  with(registerDoFuture("%dofuture%"), local = TRUE)
  stopifnot(foreach::getDoParName() == "doFuture2")
})
stopifnot(foreach::getDoParName() == "doSEQ")
