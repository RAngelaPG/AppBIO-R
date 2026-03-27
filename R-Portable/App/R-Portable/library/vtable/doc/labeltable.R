## ----eval=FALSE---------------------------------------------------------------
# labeltable(var,
#     ...,
#     out=NA,
#     count=FALSE,
#     percent=FALSE,
#     file=NA,
#     desc=NA,
#     note=NA,
#     note.align=NA,
#     anchor=NA)

## -----------------------------------------------------------------------------
#Include a single labelled variable to show how the values of that variable correspond to its value labels.
library(vtable)
library(sjlabelled)
data(efc)
labeltable(efc$e15relat)

## -----------------------------------------------------------------------------
#Include more than one variable to show, for each value of the first, what values of the others are present in the data.
data(mtcars)
labeltable(mtcars$cyl,mtcars$carb,mtcars$am)

## -----------------------------------------------------------------------------
library(vtable)
data(LifeCycleSavings)
labeltable(efc$e15relat, count = TRUE, percent = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# library(vtable)
# 
# data(LifeCycleSavings)
# labeltable(efc$e15relat,file='lifecycle_variabledocumentation')

