## ----eval=FALSE---------------------------------------------------------------
# vtable(data,
#        out=NA,
#        file=NA,
#        labels=NA,
#        class=TRUE,
#        values=TRUE,
#        missing=FALSE,
#        index=FALSE,
#        factor.limit=5,
#        char.values=FALSE,
#        data.title=NA,
#        desc=NA,
#        note=NA,
#        anchor=NA,
#        col.width=NA,
#        col.align=NA,
#        align=NA,
#        note.align='l',
#        fit.page=NA,
#        summ=NA,
#        lush=FALSE,
#        opts=list())

## -----------------------------------------------------------------------------
library(vtable)

#Example 1, using base data LifeCycleSavings
data(LifeCycleSavings)
vtable(LifeCycleSavings, out='kable')

## -----------------------------------------------------------------------------
#Example 2, using efc data with embedded variable labels
library(sjlabelled)
data(efc)
#Don't forget the handy shortcut vt()!
vt(efc)

## ----eval = FALSE-------------------------------------------------------------
# data(LifeCycleSavings)
# vtable(LifeCycleSavings)
# vtable(LifeCycleSavings,out='browser')
# vtable(LifeCycleSavings,out='viewer')
# htmlcode <- vtable(LifeCycleSavings,out='htmlreturn')
# vartable <- vtable(LifeCycleSavings,out='return')
# 
# #I can easily \input this into my LaTeX doc:
# vt(LifeCycleSavings,out='latex',file='mytable1.tex')

## ----eval=FALSE---------------------------------------------------------------
# data(LifeCycleSavings)
# vt(LifeCycleSavings,file='lifecycle_variabledocumentation')

## -----------------------------------------------------------------------------
#Note that LifeCycleSavings has five variables
data(LifeCycleSavings)
#These variable labels are taken from help(LifeCycleSavings)
labs <- c('numeric aggregate personal savings',
    'numeric % of population under 15',
    'numeric % of population over 75',
    'numeric real per-capita disposable income',
    'numeric % growth rate of dpi')
vtable(LifeCycleSavings,labels=labs)

## -----------------------------------------------------------------------------
labs <- c('numeric aggregate personal savings',NA,NA,NA,NA)
vtable(LifeCycleSavings,labels=labs)

## -----------------------------------------------------------------------------
#Note that LifeCycleSavings has five variables
#with names 'sr', 'pop15', 'pop75', 'dpi', and 'ddpi'
data(LifeCycleSavings)
#These variable labels are taken from help(LifeCycleSavings)
labs <- data.frame(nonsensename1 = c('sr', 'pop15', 'pop75'),
nonsensename2 = c('numeric aggregate personal savings',
    'numeric % of population under 15',
    'numeric % of population over 75'))
vt(LifeCycleSavings,labels=labs)

## -----------------------------------------------------------------------------
#Note that LifeCycleSavings has five variables
#with names 'sr', 'pop15', 'pop75', 'dpi', and 'ddpi'
data(LifeCycleSavings)
#These variable labels are taken from help(LifeCycleSavings)
labs <- data.frame(sr = 'numeric aggregate personal savings',
    pop15 = 'numeric % of population under 15',
    pop75 = 'numeric % of population over 75')
vtable(LifeCycleSavings,labels=labs)

## -----------------------------------------------------------------------------
data(LifeCycleSavings)
vtable(LifeCycleSavings,values=FALSE)
vtable(LifeCycleSavings)

#CO2 contains factor variables
data(CO2)
vtable(CO2)

## -----------------------------------------------------------------------------
#efc contains labeled values
#Note that the original value labels do not easily tell you what numerical
#value each label maps to, but vtable() does.
library(sjlabelled)
data(efc)
vtable(efc)

## ----eval=FALSE---------------------------------------------------------------
# data(USJudgeRatings)
# USJudgeRatings$Judge <- row.names(USJudgeRatings)
# USJudgeRatings$SecondCharacter <- 'Less Interesting'
# USJudgeRatings$ThirdCharacter <- 'Less Interesting Still!'
# 
# #Show values for all character variables
# vtable(USJudgeRatings,char.values=TRUE)
# #Or just for a subset
# vtable(USJudgeRatings,char.values=c('Judge','SecondCharacter'))

## ----eval=FALSE---------------------------------------------------------------
# library(vtable)
# 
# data(LifeCycleSavings)
# vtable(LifeCycleSavings)
# vtable(LifeCycleSavings,data.title='Intercountry Life-Cycle Savings Data',
#     desc='omit')
# vtable(LifeCycleSavings,data.title='Intercountry Life-Cycle Savings Data',
#     desc='Data on the savings ratio 1960–1970. omit')
# vtable(LifeCycleSavings,data.title='Intercountry Life-Cycle Savings Data',
#     desc='Data on the savings ratio 1960–1970',
#     note='Data from Belsley, Kuh, and Welsch (1980)')

## ----eval=FALSE---------------------------------------------------------------
# library(sjlabelled)
# data(efc)
# #The variable names in this data set are pretty short, and the value labels are
# #a little cramped, so let's move that over.
# vtable(efc,col.width=c(10,10,40,40))

## ----eval = FALSE-------------------------------------------------------------
# library(sjlabelled)
# data(efc)
# vtable(efc,col.align = 'right')

## ----eval = FALSE-------------------------------------------------------------
# library(sjlabelled)
# data(efc)
# vtable(efc,align = 'p{.3\\textwidth}cc', fit.page = '\\textwidth', out = 'latex')

## -----------------------------------------------------------------------------
library(sjlabelled)
data(efc)

vtable(efc,summ=c('mean(x)','countNA(x)'))

