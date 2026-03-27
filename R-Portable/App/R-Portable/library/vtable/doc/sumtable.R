## -----------------------------------------------------------------------------
library(vtable)
st(iris)

## ----eval=FALSE---------------------------------------------------------------
# sumtable(data,
#          vars=NA,
#          out=NA,
#          file=NA,
#          summ=NA,
#          summ.names=NA,
#          add.median=FALSE,
#          group=NA,
#          group.long=FALSE,
#          group.test=FALSE,
#          group.weights =NA,
#          col.breaks=NA,
#          digits=2,
#          fixed.digits=FALSE,
#          numformat = formatfunc(digits = digits, big.mark = ''),
#          skip.format = c('notNA(x)','propNA(x)','countNA(x)'),
#          factor.percent=TRUE,
#          factor.counts=TRUE,
#          factor.numeric=FALSE,
#          logical.numeric=FALSE,
#          logical.labels=c('No','Yes'),
#          labels=NA,
#          title='Summary Statistics',
#          note = NA,
#          anchor=NA,
#          col.width=NA,
#          col.align=NA,
#          align=NA,
#          note.align='l',
#          fit.page=NA,
#          simple.kable=FALSE,
#          obs.function=NA)
#          opts=list())

## ----eval = FALSE-------------------------------------------------------------
# data(LifeCycleSavings)
# st(LifeCycleSavings, vars = c('pop15','pop75'))

## ----eval = FALSE-------------------------------------------------------------
# data(LifeCycleSavings)
# sumtable(LifeCycleSavings)
# vartable <- vtable(LifeCycleSavings,out='return')
# 
# #I can easily \input this into my LaTeX doc:
# vt(LifeCycleSavings,out='latex',file='mytable1.tex')

## ----eval=FALSE---------------------------------------------------------------
# data(LifeCycleSavings)
# st(LifeCycleSavings,file='lifecycle_summary')

## -----------------------------------------------------------------------------
sumtable(iris,
         summ=c('notNA(x)',
                'mean(x)',
                'median(x)',
                'propNA(x)'))

## -----------------------------------------------------------------------------
#Getting complex
st(iris, col.breaks = 4,
   summ = list(
     c('notNA(x)','mean(x)','sd(x^2)','min(x)','max(x)'),
     c('notNA(x)','mean(x)')
   ),
   summ.names = list(
     c('N','Mean','SD of X^2','Min','Max'),
     c('Count','Percent')
   ))

## -----------------------------------------------------------------------------
st(iris, group = 'Species', group.test = TRUE)

## -----------------------------------------------------------------------------
st(iris, group = 'Species', group.long = TRUE)

## -----------------------------------------------------------------------------
#Let's put species in a column by itself
#There are five variables here, Species is last,
#so break the column after the first four variables.
st(iris, col.breaks = 4)

## -----------------------------------------------------------------------------
#Why not three columns?
sumtable(mtcars, col.breaks = c(4,8))

## -----------------------------------------------------------------------------
st(iris, digits = 5)

## -----------------------------------------------------------------------------
st(iris, digits = 3, fixed.digits = TRUE, numformat = NA)

## -----------------------------------------------------------------------------
st(iris, numformat = c('|cm', 'Sepal.Width' = 'percent'))

## -----------------------------------------------------------------------------
st(iris, factor.percent = FALSE, factor.counts = FALSE)

## -----------------------------------------------------------------------------
st(iris, factor.numeric = TRUE)

## -----------------------------------------------------------------------------
#Note that LifeCycleSavings has five variables
data(LifeCycleSavings)
#These variable labels are taken from help(LifeCycleSavings)
labs <- c('numeric aggregate personal savings',
    'numeric % of population under 15',
    'numeric % of population over 75',
    NA,
    'numeric % growth rate of dpi')
sumtable(LifeCycleSavings,labels=labs)

## -----------------------------------------------------------------------------
#Note that LifeCycleSavings has five variables
#with names 'sr', 'pop15', 'pop75', 'dpi', and 'ddpi'
labs <- data.frame(nonsensename1 = c('sr', 'pop15', 'pop75'),
nonsensename2 = c('numeric aggregate personal savings',
    'numeric % of population under 15',
    'numeric % of population over 75'))
st(LifeCycleSavings,labels=labs)

## -----------------------------------------------------------------------------
labs <- data.frame(sr = 'numeric aggregate personal savings',
    pop15 = 'numeric % of population under 15',
    pop75 = 'numeric % of population over 75')
sumtable(LifeCycleSavings,labels=labs)

## ----eval=FALSE---------------------------------------------------------------
# #The variable names in this data set are pretty short, and the value labels are
# #a little cramped, so let's move that over.
# st(LifeCycleSavings,
#    col.width=c(9,rep(13,7)))

## ----eval = FALSE-------------------------------------------------------------
# st(LifeCycleSavings,col.align = 'right')

## ----eval = FALSE-------------------------------------------------------------
# sumtable(iris,align = 'p{.3\\textwidth}ccccccc', fit.page = '\\textwidth', out = 'latex')

