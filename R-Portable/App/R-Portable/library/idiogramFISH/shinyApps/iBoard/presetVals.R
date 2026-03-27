#
#   objects to load
#
sorry <- "Sorry, no pandoc > 2.11, or no rmarkdown, try in Rstudio, or installing pandoc and rmarkdown"

emptydata.frame <- data.frame()
# setwd("/home/fernando/GoogleDrive/gitlab/idiogramFISH/inst/shinyApps/iBoard")
paramVec    <- readRDS("www/paramVec.rds")

iniLen      <- length(paramVec$addOTUNameVec)

paramValues <- readRDS("www/paramValues.rds")
# paramValues$alpha_valDefault <- 1
# saveRDS(paramValues, "www/paramValues.rds")
exampleVec <- c("3.1" = 1, "3.2" = 2, "3.3" = 11, "3.4" = 12, "4.2" = 13, "4.3" = 14, "4.3b" = 5, "4.5" = 15,
  "5.1" = 3, "5.2" = 4, "6.1" = 16, "6.2" = 17, "6.3" = 18, "6.4" = 19, "7.1" = 6, "7.2" = 7,
  "7.3a" = 20, "7.3b" = 21, "8.1" = 22, "8.2" = 23, "9.1" = 24, "9.2" = 8, "9.3" = 9, "9.4" = 10,
  "10.2" = 25, "10.3" = 26, "10.4" = 27, "11.1" = 28, "11.2" = 29, "11.7" = 30)

maxEx       <- max(exampleVec)

anchorHsizeFDesc     <- "`anchorHsizeF`: numeric, factor to modify horizontal size of anchor`1` (default)."
classChrNameUpDesc   <- '`classChrNameUp`: character, name of "chromosome" for col. `"chrNameUp"`. Defaults to `"Type"`'
anchorVsizeFDesc     <- "`anchorVsizeF` numeric, factor to modify vertical size of anchor `0.5` (default). Size itself is equal to `karHeiSpace`"
forbiddenMarkDesc    <- "`forbiddenMark`: character, character string or vector with mark names to be removed from plot. Not the marks but the labels."
nsmallDesc           <- "`nsmall`: numeric, rounding decimals for `chrSize` parameter. Defaults to `1`"
miniTickFactorDesc   <- "`miniTickFactor`: numeric, number of minor ticks for each labeled tick. See `useMinorTicks`. Defaults to `10`"
useMinorTicksDesc    <- "`useMinorTicks`: boolean, display minor ticks between labeled ticks in ruler. See `miniTickFactor`. Defaults to `FALSE`. (ticks without label)"
chrNameUpDesc        <- "`chrNameUp`: boolean, when `TRUE` adds secondary chromosome name from col. `chrNameUp` over chrs."

collapseCenDesc      <- "`collapseCen`: boolean, avoid spacing in ruler between short arm and long arm."
chrIdPatternRemDesc  <- "`chrIdPatternRem`: character, regex pattern to remove from chr. names"
classChrNameDesc     <- '`classChrName`: character, name of "chromosome" when in micrometers (apparently). Defaults to `"Chr."`.'
groupSeparDesc       <- "`groupSepar`: numeric, factor for affecting chr. spacing `chrSpacing` among groups"
classGroupNameDesc   <- '`classGroupName`: character, name of groups. Defaults to `""`'
bMarkNameAsideDesc   <- '`bMarkNameAside`: boolean, when `TRUE` and `legend="inline"`, shows marks in `bannedMarkName` as `legend="aside"` would do.'
markNewLineDesc      <- "`markNewLine`: character, character to split mark Names into different lines. Applies to `square` marks."
dfChrSizeDesc        <- "`dfChrSize`: mandatory data.frame or .csv file name, with columns: `OTU` (optional), `chrName` (mandatory), `shortArmSize`, `longArmSize` for monocen. or `chrSize` for holocen."
dfMarkPosDesc        <- "`dfMarkPos`: name of the data.frame of positions of marks or .csv file. Includes GISH and centromeric marks (cen). columns: `OTU` (opt), `chrName`, `markName` (name of site), `chrRegion` (for monocen. and opt for whole arm (`w`) in holocen.), `markDistCen` (for monocen.), `markPos` (for holocen.), `markSize`; column `chrRegion`: use `p` for short arm, `q` for long arm, `cen` for centromeric mark and `w` for whole chr. mark; column `markDistCen`: use distance from centromere to mark, not necessary for cen. marks (`cen`), `w`, `p`, `q` (when whole arm). See also param. `markDistType`"
dfMarkColorDesc      <- '`dfMarkColor`: name of the data.frame or .csv file of marks characteristics. Optional. Specifying colors and style for marks (sites); columns: `markName`, `markColor`, `style`. style accepts: `"square"`,`"squareLeft"`, `"dots"`, `"cM"`, `"cMLeft"`,`"cenStyle"`, `"upArrow"`, `"downArrow"`, `"exProtein"`. (if column `style` missing all (except `5S`) are plotted as in param. `defaultStyleMark`).'

notesDesc            <- "`notes`: data.frame with columns `OTU` and `note` for adding notes to each OTU, they appear to the right of the karyotype"
leftNotesDesc <- "`leftNotes`: data.frame with columns `OTU` and `note` for adding notes to each OTU, they appear to the left of the karyotype"
leftNotesUpDesc <- "`leftNotesUp`: data.frame, (to the left), similar to `leftNotes`, but intended for placement over kar."
originDesc           <- '`origin`: (`"b"`) If you measure your mark from the bottom of chromosome use `origin = "b"`, or `"t"` from top. Applies to holocentrics. (monocentrics marks are measured from centromere)'
OTUfamilyDesc        <- "`OTUfamily`: character, font family for OTU name."
xModMonoHoloRateDesc <- "`xModMonoHoloRate`:	numeric, factor to shrink `xModifier` for holocen. 5 means 5 times smaller (quotient)."
chrLabelSpacingDesc  <- "`chrLabelSpacing`: numeric, for `circularPlot=TRUE`. Spacing of chr. names. Defaults to `0.5`"
labelSpacingDesc <- "`labelSpacing`: numeric, for `circularPlot=TRUE`. Spacing of mark labels. Defaults to `0.7`"
rotationDesc <- "`rotation`: numeric, anti-clockwise rotation, defaults to `0.5` which rotates first chr. from top to -90 degrees. (-0.5*π = 9 o'clock)"
labelOutwardsDesc <- "`labelOutwards`: boolean, inline labels projected outwards"
shrinkFactorDesc <- "`shrinkFactor`:	numeric, for `circularPlot=TRUE` percentage of usage of circle. Defaults to `0.9`"
radiusDesc <- "`radius`: numeric, for `circularPlot=TRUE`. Affects radius of karyotypes. Defaults to `0.5`"
separFactorDesc <- "`separFactor`: numeric, for `circularPlot=TRUE` modify separation of concentric karyotypes. Defaults to `1.5`"
circleCenterDesc <- "`circleCenter`:	numeric, for `circularPlot=TRUE`. Affects coordinates of center of circles. Affects `legend='aside'` position."
OTUplacingDesc <- '`OTUplacing`:	character, for `circularPlot=TRUE`. location of OTU name. Defaults to `"first"`, which plots name of OTU near first chr. `"number"` places number near 1st chr. and index and name of OTU to the right or center. `"simple"` places name of OTU to the right or center without numbering. See also `OTUcentered`'
OTUsrtDesc <- "`OTUsrt`: numeric, for `circularPlot=TRUE` and `OTUplacing='first'` Angle to use for OTU names. Defaults to `0`. See `OTUplacing`"
OTUjustifDesc <- '`OTUjustif`: numeric, for `circularPlot=TRUE` and `OTUplacing="number"` or `"simple"`. Justification of OTU name. `0` = left (Default); use `0.5` for centered. See `?text` -> `adj`'
OTULabelSpacerxDesc <- '`OTULabelSpacerx`: numeric, for `circularPlot=TRUE` and `OTUplacing="number"` or `"simple"`. Modifies x names position'
OTUlegendHeightDesc <- '`OTUlegendHeight`: numeric, for `circularPlot=TRUE` and `OTUplacing="number"` or `"simple"`. Modifies y names separation'
legendDesc <- '`legend`: (`"aside"`) If you wanto to plot the names of marks near each chromosome use `legend = "inline"`, to the right of karyotypes use `legend = "aside"`, otherwise use `legend = ""` for no legend. See `markLabelSpacer`'
remSimiMarkLegDesc <- '`remSimiMarkLeg`: boolean, when `legend="aside"`, if you use the `pattern` arg., you can get several marks with "same" name. When `TRUE` this remove this pseudoduplicates from legend. Be sure that this pseudoduplicates have the same color, otherwise you should use `FALSE` (default).'
legendWidthDesc <- '`legendWidth`: (`1.7`) numeric, factor to modify the width of the square and dots of legend. For `legend="aside"`.'
patternDesc <- '`pattern`: (`""`) REGEX pattern to eliminate from the marks name when plotting. See human karyotype chapter for example.'
legendHeightDesc <- '`legendHeight`: (`NA`) numeric, factor to modify the height of the square and dots of legend. For `legend="aside"`.'
bannedMarkNameDesc <- "`bannedMarkName`: character, character string or vector with mark names to be removed from plot. Not the marks but the labels. See `bMarkNameAside`"
markLabelSizeDesc <- "`markLabelSize`: (`1`) Determines the size of text of the legend."
markLabelSpacerDesc <- "`markLabelSpacer`: (`1`) When `legend = \"aside\"` determines the separation of legends from the karyotype right side"
legendYcoordDesc <- '`legendYcoord`: numeric, modify Y position of legend when `legend="aside"`'
OTUasNoteDesc <- "`OTUasNote`: (`FALSE`) See also `notes`. If `TRUE` OTU name is written to the right, as `notes`."
parseStr2langDesc <- '`parseStr2lang`: boolean, parse string in `notes` with function `str2lang(paste0("paste(",note,")") )` for ex: `"italic(\'C. sinensis\'), \' Author\'"`. See `notes`, `leftNotes`,`leftNotesUp`. '
notesTextSizeDesc <- "`notesTextSize`: (`0.4`) numeric, font size of notes, see `notes`"
noteFontDesc <- "`noteFont`: numeric  `1` for normal,  `2` for bold,  `3` for italics,  `4` for bold-italics. See `notes`"
notesPosXDesc <- "`notesPosX`: (`0.5`) numeric, moves right notes in the x axis"
notesPosYDesc <- "`notesPosY`: (`0.5`) numeric, moves right notes in the y axis"
leftNotesTextSizeDesc <- "`leftNotesTextSize`: (`0.4`) numeric, font size of notes, see `leftNotes`"
leftNoteFontDesc <- "`leftNoteFont`: numeric  `1` for normal,  `2` for bold,  `3` for italics,  `4` for bold-italics. See `leftNotes`"
leftNotesPosXDesc <- "`leftNotesPosX`: (`0.5`) numeric, moves left notes in the x axis"
leftNotesPosYDesc <- "`leftNotesPosY`: numeric, move `leftNotes` down or up (y axis)"
leftNotesUpTextSizeDesc <- "`leftNotesUpTextSize`: (`0.4`) numeric, font size of notes, see `leftNotesUp`"
leftNoteFontUpDesc <- "`leftNoteFontUp`: numeric  `1` for normal,  `2` for bold,  `3` for italics,  `4` for bold-italics. See `leftNotesUp`"
leftNotesUpPosXDesc <- "`leftNotesUpPosX`: numeric, move `leftNotesUp`  (x axis)"
leftNotesUpPosYDesc <- "`leftNotesUpPosY`: numeric, move `leftNotesUp` down or up (y axis)"
protrudingDesc <- "`protruding`: numeric, when style of mark is `cM`, fraction of chrWidth to stretch marker. Defaults to `0.2`"
arrowheadDesc <- "`arrowhead`: numeric, proportion of head of arrow (mark styles: `upArrow`,`downArrow`). Defaults to `0.3`"
useOneDotDesc <- "`useOneDot`: boolean, use one dot instead of two in style of marks dots. Defaults to `FALSE`"
cMBeginCenterDesc <- "`cMBeginCenter`: boolean, start position of `cM` and `cMLeft` marks. If `TRUE`, starts in the center (width) of chr. . Defaults to `FALSE`"
pMarkFacDesc <- "`pMarkFac` numeric, fraction of chr. size for `exProtein` style marks. Defaults to `0.25`"
markDistTypeDesc <- '`markDistType`: (`"beg"`) If you measure your marks to the beginning of mark use  `markDistType = "beg"`, if to the center of the mark, use `"cen"`. '
hideCenLinesDesc <- "`hideCenLines`: numeric, factor to multiply line width (lwd) used for covering cen. border, when `chrColor` is `white` or when `gishCenBorder=TRUE`"
lwd.marksDesc <- "`lwd.marks`: thickness of most marks. Except `cM` marks and centr. related marks. See `lwd.chr`, `lwd.cM`. Defaults to `lwd.chr` when `99`"
lwd.mimicCenDesc <- "`lwd.mimicCen`: thickness of lines of `cenStyle` marks; affects only lateral borders. Defaults to `lwd.chr`"
lwd.cMDesc <- "`lwd.cM`: thickness of cM marks. Defaults to `lwd.chr`"
addMissingOTUAfterDesc <- "`addMissingOTUAfter`: (`NA`) character vector, Pass to this parameter a vector of OTUs after which empty spaces (ghost karyotypes) must be added. See `missOUTspacings` and the `phylogeny` chapter"
missOTUspacingsDesc <- "`missOTUspacings`: (`0`) numeric vector. With the same length of `addMissingOTUAfter`. Pass to this parameter the number of ghost karyotypes following the OTUs passed to `addMissingOTUAfter`. See the `phylogeny` chapter"
specialOTUNamesDesc <- "`specialOTUNames`: character vector, normally title of ruler is `μm` or `Mb` (big numbers). Use this param. to be able to put a different unit in ruler title. See `specialyTitle`"
specialyTitleDesc <- "`specialyTitle`: character, title of ruler if OTU is in `specialOTUNames`. Will not apply if `MbThreshold` met. In that case use `MbUnit`"
specialChrWidthDesc <- "`specialChrWidth`: numeric, relative chromosome width. Defaults to `0.5` for OTUs in `specialOTUNames`"
specialChrSpacingDesc <- "`specialChrSpacing` numeric, horizontal spacing among chromosomes for OTUs in `specialOTUNames`, see also  `chrWidth`. Defaults to `0.5` "
addOTUNameDesc <- "`addOTUName`: (`TRUE`) If `TRUE` adds name of species (OTU) under karyotype"
OTUTextSizeDesc <- "`OTUTextSize`: (`1`) font size of OTU names, except when `OTUasNote=TRUE` see `notesTextSize`"
OTUfontDesc <- "`OTUfont`: numeric, `1` for normal, `2` for bold, `3` for italics, `4` for bold-italics"
circularPlotDesc <- "`circularPlot`:	boolean, if `TRUE` chromosomes are plotted in concentric circles. Defaults to `FALSE`. See `verticalPlot`"
orderChrDesc <- '`orderChr`: (`size`) character, when `"size"`, sorts chromosomes by total length from the largest to the smallest. `"original"`: preserves d.f. order. `"name"`: sorts alphabetically; `"group"`: sorts by group name; `"chrNameUp"`: sorts according to column `chrNameUp`. See `chrNameUp`'
chrIdDesc <- '`chrId`: (`"original"`) If you want to rename chromosomes from 1 to n use `chrId = "simple"`. For original names use `chrId = "original"`. For no names use  `chrId = ""`'
chrWidthDesc <- "`chrWidth`: (`0.5`) Determines the width of chromosomes"
chrSpacingDesc <- "`chrSpacing`: (`0.5`) Determines the horizontal spacing among chromosomes"
lwd.chrDesc <- "`lwd.chr`: (`0.5`) width of border lines for chr. and marks when related param. absent."
karSeparDesc <- "`karSepar`: (`TRUE`) If `TRUE` reduces the space among karyotypes. `FALSE` = equally sized karyotypes or `TRUE` = equally spaced karyotypes. Incompatible with `addMissingOTUAfter`"
verticalPlotDesc <- "`verticalPlot`: boolean, when `TRUE` karyotypes are plotted vertically, otherwise, horizontally. Defaults to `TRUE`"
karHeightDesc <- "`karHeight`: (`2`) Vertical size of karyotypes considering only chromosomes. for ex `karHeight = 1`"
karHeiSpaceDesc <- "`karHeiSpace`: (`2.5`) Vertical size of karyotypes including spacer. for ex `karHeiSpace = 1.2`. Use with `karSepar=FALSE`"
amoSeparDesc <- "`amoSepar`: (`9`) For `karSepar = TRUE`, if zero, no space among karyotypes. Amount of separation.  if overlap, increase this and `karHeiSpace`"
karSpaceHorDesc <- "`karSpaceHor`: numeric, separation among horizontal karyotypes. When `verticalPlot=FALSE`. Defaults to `0`"
chromatidsDesc <- "`chromatids`: boolean, when `TRUE` shows separated chromatids. Defaults to `TRUE`"
holocenNotAsChromatidsDesc <- "`holocenNotAsChromatids`:	boolean, when `TRUE` and `chromatids=TRUE` does not plot holocen kar. with chromatids. Defaults to `FALSE`."
xModifierDesc <- "`xModifier`: numeric, for `chromatids=TRUE`, separation among chromatids. Quotient for `chrWidth`. Defaults to `12 = chrWidth/12`"
rulerDesc <- "`ruler`: (`TRUE`) When `TRUE` displays ruler to the left of karyotype, when `FALSE` shows no ruler"
ceilingFactorDesc <- "`ceilingFactor`: (`0`) numeric, affects number of decimals for ceiling. Affects max. value of ruler. When `threshold` is greater than `35` this may have to be negative. "
rulerPosDesc <- '`rulerPos`: (`-0.5`) Absolute position of ruler, corresponds to "pos" argument of the function `axis` of R plots'
xPosRulerTitleDesc <- "`xPosRulerTitle`: (`2.6`) Modifies the horizontal position of the title of rulers (Mb, etc). Moves to left from 1st chr. in `chrSpacing` times"
rulerNumberPosDesc <- "`rulerNumberPos`: (`0.5`) numeric, Modify position of numbers of ruler"
yTitleDesc <- "`yTitle`: (`μm`) character, units for common title."
rulerIntervalDesc <- "`rulerInterval`: numeric, intervals in ruler."
rulerIntervalcMDesc <- "`rulerIntervalcM`: numeric, intervals in ruler of OTU in `specialOTUNames`."
rulerIntervalMbDesc <- "`rulerIntervalMb`: numeric, intervals in ruler of OTU with data in Mb (>`MbThreshold`) and absent from `specialOTUNames`."
ruler.tckDesc <- '`ruler.tck`: (`-0.02`) tick size of ruler, corresponds to "tck" argument of `axis` function'
thresholdDesc <- "`threshold`: (`35`) This is the max. value allowed for the main two significative digits, otherwise scale will shrink. For example, after 35 μm (Default), apparent size will be 3.5 (not 35) and scale interval will change. See `ceilingFactor`: you may have to use `-1` for it. Introduced in 1.13"
rulerNumberSizeDesc <- "`rulerNumberSize`: (`1`) Size of number's font in ruler"
rulerTitleSizeDesc <- "`rulerTitleSize`: numeric font size of units of ruler."
autoCenSizeDesc <- "`autoCenSize`: boolean, when `TRUE` ignores `centromereSize`. Affected heavily by `threshold`"
cenFormatDesc <- '`cenFormat`: boolean, when "triangle", cen. has triangular aspect. When "rounded", it has rounded aspect (Default). "inProtein" for using the mark with style of same name.'
cenFactorDesc <- "`cenFactor`: numeric, modifies any cen. mark and cen. size. Defaults to `1`"
centromereSizeDesc <- "`centromereSize`: Apparent size of centromeres. Requires `autoCenSize = FALSE`"
squarenessDesc <- '`squareness`: (`4`) Squared or rounded vertices when marks of the "square" style (defined in data.frame passed to `dfMarkColor`). Affects chromosomes also. Smaller numbers = more rounded'
markNDesc <- "`markN`: numeric vertices number for round corners of marks"
nDesc <- "`n`: numeric vertices number for round corners of chr."
ylimBotModDesc <- "`ylimBotMod`: (`0.2`) modify `ylim` bottom component of plot adding more space"
ylimTopModDesc <- "`ylimTopMod`: (`0.2`) modify `ylim` top component of plot adding more space."
xlimLeftModDesc <- '`xlimLeftMod`: (`1`) modifies `xlim` left (first) component of the plot as in any "R-plot"'
xlimRightModDesc <- "`xlimRightMod`: (`2`) `xlim` (right) modification by adding space to the right of idiograms"
moveAllKarValueHorDesc <- "`moveAllKarValueHor`: numeric, similar to `mkhValue`, but affects all karyotypes."
moveAllKarValueYDesc <- "`moveAllKarValueY`: numeric, similar to `moveAllKarValueHor`, but affects y axis."
indexIdTextSizeDesc <- "`indexIdTextSize`: numeric, font size of chr. and kar. indices and chromosome name. Defaults to `1`"
distTextChrDesc <- "`distTextChr`: Vertical distance from indices (text) to the chromosome."
morphoDesc <- '`morpho`: (`"both"`) character, if `"both"` (default) prints the Guerra (1986) and Levan (1964) classif. of cen. position.  , use also `"Guerra"` or `"Levan"` or `""` for none. See `?armRatioCI` also (function).'
chrIndexDesc <- '`chrIndex`: (`"both"`) character, add arm ratio with `"AR"` and centromeric index with `"CI"`, or `"both"` (Default), or `""` for none to each chromosome [@Levan1964]. See `armRatioCI`also.'
chrSizeDesc <- "`chrSize`: boolean, when `TRUE` adds total chr size under each chr. Defaults to `FALSE`"
karIndexDesc <- "`karIndex`: (`TRUE`) boolean. Adds karyotype indices A (intra - cen) and A2 (inter - size) [@Watanabe1999; @Zarco1986new]. Disable with `karIndex = FALSE`"
chrSizeMbpDesc <- "`chrSizeMbp`: boolean, when `TRUE` adds total Mbp chr. size to each chr. provided, there is a `Mbp` column in `dfChrSize` data.frame. Defaults to `FALSE`. If data in columns `shortArmSize`, or col. `chrSize` is in millions ('Mbp'). Use `chrSize=TRUE` not this one (not column `Mbp`, you don't need this)."
karIndexPosDesc <- "`karIndexPos`: (`0.5`) numeric. Move karyotype index."
nameChrIndexPosDesc <- "`nameChrIndexPos`: (`2`) Modify position of name (CI, AR) of chr. indices, numeric value."
perAsFractionDesc <- "`perAsFraction`: boolean, when `TRUE` % is shown as fraction. Defaults to `FALSE`. See `markPer`"
markPerDesc <- "`markPer`: character, name of marks to calculate % of mark in chr. and add it to plot. See `perAsFraction`"
showMarkPosDesc <- "`showMarkPos`: boolean, adds position of marks under karyotype (fraction 0-1) when `TRUE`. Defaults to `FALSE`"
bToRemoveDesc <- "`bToRemove`: character, bands to remove from calc. of pos, when `showMarkPos = TRUE`"
anchorDesc <- "`anchor`: boolean, when TRUE, plots a parent progeny structure in karyotypes in `moveKarHor`"
moveKarHorDesc <- "`moveKarHor`: character, OTUs' names of karyotypes that should be moved horizontally. See `mkhValue`"
mkhValueDesc <- "`mkhValue`: numeric, value to move kar. hor. See `moveKarHor`"
karAnchorLeftDesc <- "`karAnchorLeft`: character, OTUs' names of karyotypes to the right of your desired anchor. For `verticalPlot=FALSE`"
moveAnchorVDesc <- "`moveAnchorV`: numeric, displace anchor vertical portion to right or left. See `anchor`"
moveAnchorHDesc <- "`moveAnchorH`: numeric, displace anchor horizontal portion to right or left. See `anchor`"
cenColorDesc <- "`cenColor`: Determines the color of centromeres. if GISH use `NULL`. Defaults to `chrColor`"
chrColorDesc <- '`chrColor`: (`"gray"`) Determines the color of chromosomes'
colorBorderMarkDesc <- "`colorBorderMark`: character, without default, pass a name of a color to use as border of marks. See `borderOfWhiteMarks`"
chrBorderColorDesc <- "`chrBorderColor`: character, color for border of chromosomes, defaults to `chrColor`"
fixCenBorderDesc <- "`fixCenBorder`: boolean, when `TRUE` uses `chrColor` as centromere (and cen. mark) border color. See also `cenColor`, `chrColor`, `colorBorderMark`, `borderOfWhiteMarks`. No default value."
gishCenBorderDesc <- "`gishCenBorder`: boolean, when `TRUE`, cen. mark border color is the same as mark color, ignoring `colorBorderMark`. No default."
mycolorsDesc <- "`mycolors`: optional, character vector with colors' names, which are associated automatically with marks according to their order in the data.frame of position of marks. See this ordering with `unique(dfMarkPos$markName)`. Argument example: `mycolors = c(\"red\",\"chartreuse3\",\"dodgerblue\")`. Not mandatory for plotting marks, package has default colors."
alpha_valDesc <- "`alpha_val`: alpha_val numeric vector, make marks transparent, accepts values from 0 to 1, see `scales::alpha`"

helpString <- "Start in pages (left): Examples, Nucleotides or data.frames
              \n modify browser zoom with Ctrl [+/-]
              \n If you are in Rstudio Desktop Viewer, click its button 'Open Browser'
              for better memory management"

presetUseText <- "The preset file contains all, use it in the 'Examples' page"
