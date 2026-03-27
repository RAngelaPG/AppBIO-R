## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
options(rmarkdown.html_vignette.check_title = FALSE)
library(statgenHTP)
library(ggplot2)

## ----importPhenovator---------------------------------------------------------
data("PhenovatorDat1")

## ----headPhenovato, echo=FALSE, message=FALSE---------------------------------
knitr::kable(head(PhenovatorDat1), align=c('c','c'), booktabs = TRUE)

## ----createTP-----------------------------------------------------------------
## Create a TP object containing the data from the Phenovator.
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            experimentName = "Phenovator",
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            repId = "Replicate",
                            plotId = "pos",
                            rowNum = "y", colNum = "x",
                            addCheck = TRUE,
                            checkGenotypes = c("check1", "check2", "check3", "check4"))
summary(phenoTP)

## ----getTimepoints------------------------------------------------------------
## Extract the time points table.
timepoint <- getTimePoints(phenoTP)

## ----getTimepointsbis, echo=FALSE, message=FALSE------------------------------
knitr::kable(head(timepoint), align=c('c','c'), padding = 0)

## ----layoutPlot, fig.height=4, fig.width=5, fig.align = 'center'--------------
## Plot the layout for the third time point.
plot(phenoTP, 
     plotType = "layout",
     timePoints = 3)

## ----layoutPlotHL, fig.height=4, fig.width=5, fig.align = 'center'------------
## Plot the layout for the third time point with the check genotypes highlighted.
plot(phenoTP, 
     plotType = "layout",
     timePoints = 3,  
     highlight = c("check1", "check2", "check3", "check4"))

## ----layoutPlotSG, fig.height=7, fig.width=8, fig.align = 'center'------------
## Plot the layout for the third time point.
plot(phenoTP, 
     plotType = "layout",
     timePoints = 3,  
     highlight = c("check1", "check2", "check3", "check4"),
     showGeno = TRUE)

## ----layoutPlotheatmap, fig.height=7, fig.width=8, fig.align = 'center'-------
## Plot the layout for the third time point.
plot(phenoTP, 
     plotType = "layout",
     timePoints = 3,  
     traits = "EffpsII")

## ----rawVator, include=TRUE, fig.height=3, fig.width=7, fig.align = 'center'----
## Create the raw data time courses for three genotypes.
plot(phenoTP, 
     traits = "EffpsII",
     plotType = "raw",
     genotypes = c("G001", "G002", "check1"))

## ----boxPlot, fig.height=4.5, fig.width=7, fig.align = 'center'---------------
## Create a boxplot for "EffpsII" using the default all time points.
plot(phenoTP, 
     plotType = "box",
     traits = "EffpsII") 

## ----corPlot, fig.height=5, fig.width=6, fig.align = 'center'-----------------
## Create a correlation plot for "EffpsII" for a selection of time points.
plot(phenoTP, 
     plotType = "cor",
     traits = "EffpsII",
     timePoints = seq(from = 1, to = 73, by = 5))

## ----setParamVator------------------------------------------------------------
# First select a subset of plants, for example here 4 plants.
plantSel <- c("c1r17","c13r17","c6r51","c21r24") 
# Then run on the subset
resuVatorHTP <- detectSingleOut(TP = phenoTP,
                                trait = "EffpsII",
                                plotIds = plantSel,
                                confIntSize = 3,
                                nnLocfit = 0.1)

## ----headOutVator, echo=FALSE, message=FALSE----------------------------------
knitr::kable(head(resuVatorHTP), align = c('c','c'), padding = 0)

## ----plotOutVatorFA1, fig.height=1.5, fig.width=5, fig.align = 'center'-------
plot(resuVatorHTP,
     outOnly = FALSE)

## ----rmSingleOutVator---------------------------------------------------------
phenoTPOut <- removeSingleOut(phenoTP,
                              resuVatorHTP)

## ----fitSp, message=FALSE-----------------------------------------------------
## Fit a model for a few time points.
modPhenoSp <- fitModels(TP = phenoTPOut, 
                        trait = "EffpsII",
                        timePoints = c(1, 33, 36, 54, 73))
summary(modPhenoSp)

## ----plotSpatPerc,  fig.height=4, fig.width=7.3, message=FALSE----------------
plot(modPhenoSp,
     timePoints = 36,
     plotType = "spatial",
     spaTrend = "percentage")

## ----plotTimeLapse,  fig.height=5, fig.width=6, message=FALSE, eval=FALSE-----
# plot(modPhenoSp,
#      plotType = "timeLapse",
#      outFile = "TimeLapse_modPhenoSp.gif")

## ----plotSpPred, message=FALSE, fig.height=2, fig.width=6, fig.align = 'center'----
plot(modPhenoSp, 
     plotType = "rawPred",
     genotypes = c("G007", "G058")) 

## ----plotSpCorr, message=FALSE, fig.height=2, fig.width=6, fig.align = 'center'----
plot(modPhenoSp, 
     plotType = "corrPred",
     genotypes = c("G007", "G058"))

## ----plotSpHerit, message=FALSE, fig.height=2.5, fig.width=3.5, fig.align = 'center'----
plot(modPhenoSp, 
     plotType = "herit",
     yLim = c(0.5, 1))

## ----plotSpVar, fig.height=3, fig.width=5, message=FALSE, fig.align = 'center'----
plot(modPhenoSp, 
     plotType = "variance")

## ----plotSpED, message=FALSE, fig.height=3, fig.width=5, fig.align = 'center'----
plot(modPhenoSp, 
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"),
     EDType = "ratio")

## ----getFun, message=FALSE----------------------------------------------------
## Extract the genotypic predictions for one time point: 
genoPredSp <- getGenoPred(modPhenoSp)

## ----getPred, echo=FALSE, message=FALSE---------------------------------------
knitr::kable(head(genoPredSp$genoPred), align = "c", padding = 0)

## ----fitSplineVator, message=FALSE, warning=FALSE-----------------------------
data(spatCorrectedVator)  
# Fit P-splines using on a subset of genotypes.
subGenoVator <- c("G160", "G151")

fit.spline <- fitSpline(inDat = spatCorrectedVator,
                        trait = "EffpsII_corr",
                        genotypes = subGenoVator,
                        knots = 50,
                        useTimeNumber = TRUE,
                        timeNumber = "timeNumHour")
# Extracting the tables of predicted values and P-spline coefficients
predDat <- fit.spline$predDat
coefDat <- fit.spline$coefDat

## ----plotSplGeno,  fig.height=4, fig.width=7, message=FALSE-------------------
plot(fit.spline,
     genotypes = "G160")

## ----OutVator, message=FALSE, warning=FALSE-----------------------------------
outVator <- detectSerieOut(corrDat = spatCorrectedVator,
                           predDat = predDat,
                           coefDat = coefDat,
                           trait = "EffpsII_corr",
                           genotypes = subGenoVator,
                           thrCor = 0.9,
                           thrPca = 30)

## ----headOutPoint, echo=FALSE, message=FALSE----------------------------------
knitr::kable(head(outVator), align = "c", booktabs = TRUE, row.names = FALSE)

## ----plotOutVator,  fig.height=6, fig.width=6, message=FALSE, warning=FALSE----
plot(outVator, genotypes = "G151")

## ----rmSerieOutVator----------------------------------------------------------
fit.splineOut <- removeSerieOut(fitSpline = fit.spline,
                                serieOut = outVator)

## ----summaryData, message=FALSE, warning=FALSE, eval=TRUE, fig.align='center'----
data(spatCorrectedVator)
spatCorrectedVator[["pop"]] <- as.factor(rep("Pop1", nrow(spatCorrectedVator)))
str(droplevels(spatCorrectedVator[spatCorrectedVator$genotype %in% subGenoVator,]))

## ----fitPsHDMVator, message=FALSE, warning=FALSE------------------------------
## Fit P-spline HDM.
fit.psHDM  <- fitSplineHDM(inDat = spatCorrectedVator,
                           genotypes = subGenoVator,
                           trait = "EffpsII_corr",
                           useTimeNumber = TRUE,
                           timeNumber = "timeNumHour",                           
                           pop = "pop",
                           genotype = "genotype",
                           plotId = "plotId",
                           weights = "wt",
                           difVar = list(geno = FALSE, plot = FALSE),
                           smoothPop = list(nseg = 20, bdeg = 3, pord = 2),
                           smoothGeno = list(nseg = 20, bdeg = 3, pord = 2),
                           smoothPlot = list(nseg = 20, bdeg = 3, pord = 2),
                           trace = FALSE)

## ----predictSplineHDMVatorNum, message=FALSE, warning=FALSE-------------------
## Predict P-spline HDM.
pred.psHDM <- predict(object = fit.psHDM,
                      newtimes = seq(min(fit.psHDM$time[["timeNumber"]]),
                                     max(fit.psHDM$time[["timeNumber"]]),
                                     length.out = 100),
                      pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                      se = list(pop = TRUE, geno = TRUE, plot = FALSE),
                      trace = FALSE)

## ----plotPredPopVator, fig.height=4, fig.width=5, message=FALSE, warning=FALSE, fig.align='center'----
plot(pred.psHDM, plotType = "popTra", themeSizeHDM = 10)

## ----plotPredGenoTraVator, fig.height=4, fig.width=5, message=FALSE, warning=FALSE, fig.align='center'----
plot(pred.psHDM, plotType = "popGenoTra", themeSizeHDM = 10)

## ----plotPredGenoDevVator, fig.height=4, fig.width=5, message=FALSE, warning=FALSE, fig.align='center'----
plot(pred.psHDM, plotType = "genoDev", themeSizeHDM = 10)

## ----plotPredPlotVator, fig.height=4, fig.width=6, message=FALSE, warning=FALSE, fig.align='center'----
plot(pred.psHDM, 
     plotType = "genoPlotTra", 
     themeSizeHDM = 10)

## ----paramVator, fig.height=2, fig.width=4, message=FALSE, warning=FALSE, fig.align='center'----
paramVator1 <- estimateSplineParameters(x = fit.splineOut,
                                        estimate = "predictions",
                                        what = "AUC",
                                        timeMin = 330,
                                        timeMax = 432,
                                        genotypes = subGenoVator)

plot(paramVator1, plotType = "box")

## ----paramPhenoGenoPred, fig.height=2, fig.width=4, message=FALSE, warning=FALSE, fig.align='center'----
paramVator2 <-
  estimateSplineParameters(x = pred.psHDM,
                           what = "min",
                           fitLevel = "plot",
                           estimate = "predictions")

plot(paramVator2, plotType = "box")

