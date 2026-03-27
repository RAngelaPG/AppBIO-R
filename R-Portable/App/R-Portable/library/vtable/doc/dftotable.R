## ----eval = FALSE-------------------------------------------------------------
# dftoHTML(data,
#          out=NA,
#          file=NA,
#          note=NA,
#          anchor=NA,
#          col.width=NA,
#          col.align=NA,
#          row.names=FALSE,
#          no.escape=NA)

## ----eval = FALSE-------------------------------------------------------------
# library(vtable)
# 
# data(LifeCycleSavings)
# dftoHTML(LifeCycleSavings)

## ----eval=FALSE---------------------------------------------------------------
# library(vtable)
# 
# data(LifeCycleSavings)
# dftoHTML(LifeCycleSavings)
# dftoHTML(LifeCycleSavings,out='browser')
# dftoHTML(LifeCycleSavings,out='viewer')
# htmlcode <- dftoHTML(LifeCycleSavings,out='htmlreturn')

## ----eval=FALSE---------------------------------------------------------------
# data(LifeCycleSavings)
# dftoHTML(LifeCycleSavings,file='lifecycledata_htmlversion.html')

## ----eval = FALSE-------------------------------------------------------------
# dftoHTML(LifeCycleSavings,note='Data from Belsley, Kuh, and Welsch 1980').

## ----eval = FALSE-------------------------------------------------------------
# #Let's make sr much bigger for some reason
# dftoHTML(LifeCycleSavings,col.width=c(60,10,10,10,10))

## ----eval=FALSE---------------------------------------------------------------
# dftoHTML(LifeCycleSavings,row.names=TRUE)

## ----eval=FALSE---------------------------------------------------------------
# #Don't escape columns 1 or 2
# dftoHTML(LifeCycleSavings,no.escape=1:2)

## ----eval = FALSE-------------------------------------------------------------
# dftoLaTeX(data,
#          file=NA,
#          frag=TRUE,
#          title=NA,
#          note=NA,
#          anchor=NA,
#          align=NA,
#          row.names=FALSE,
#          no.escape=NA)

## ----eval = FALSE-------------------------------------------------------------
# library(vtable)
# 
# data(LifeCycleSavings)
# dftoLaTeX(LifeCycleSavings)

## ----eval=FALSE---------------------------------------------------------------
# data(LifeCycleSavings)
# dftoLaTeX(LifeCycleSavings,file='lifecycledata_latexversion.tex')

## ----eval = FALSE-------------------------------------------------------------
# dftoLaTeX(LifeCycleSavings,note='Data from Belsley, Kuh, and Welsch 1980').

## ----eval = FALSE-------------------------------------------------------------
# dftoLaTeX(LifeCycleSavings,anchor='tab:LCS')

## ----eval = FALSE-------------------------------------------------------------
# dftoLaTeX(LifeCycleSavings,row.names=TRUE,align='p{.25\\textwidth}ccccc')

## ----eval=FALSE---------------------------------------------------------------
# dftoLaTeX(LifeCycleSavings,row.names=TRUE)

## ----eval=FALSE---------------------------------------------------------------
# #Don't escape columns 1 or 2
# dftoLaTeX(LifeCycleSavings,no.escape=1:2)

