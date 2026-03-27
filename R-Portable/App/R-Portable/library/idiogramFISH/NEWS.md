
<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# idiogramFISH 2.0.13

22-08-2023

Docs:

- Fix arguments descriptions

# idiogramFISH 2.0.12

06-06-2023

Docs:

- Fix link on DESCRIPTION file

# idiogramFISH 2.0.11

05-04-2023

Bug:

- Fix genBank processing when only one match

# idiogramFISH 2.0.10

05-04-2023

param:

- `alpha_val`: modify transparency of marks

# idiogramFISH 2.0.9

14-09-2022

- refactoring of code

# idiogramFISH 2.0.8

27-12-2021

- Accepts .csv files to notes parameters
- `circularPlot=TRUE` now compatible with `chromatids=TRUE` showing
  chromatids in circ. plots
- centromere is now subtracted from arm, so ruler is continuous in
  monocen. (`collapseCen`).
- new `chrRegion` possible (column in data.frame `dfMarkPos`): `pcen`
  `qcen`. It’s behavior is similar to `cen`

param:

- `dfCenMarks`: deprecated
- `collapseCen`: boolean, avoid spacing in ruler between short arm and
  long arm.
- `anchorHsizeF`: numeric, factor to modify horizontal size of anchor
  `1` (default).

Shiny:

- Change in default values for new data.frames
- Fixed parsing of `yTitle` in code tab
- Download and upload custom presets

# idiogramFISH 2.0.6

27-09-2021

- Added possibility of showing several % (span) of marks, ex.
  `markPer = c("5S","45S")`
- `perMark` and `posCalc` functions produce `data.frames`
- Added parameters: `autoCenSize`, `leftNotesUpPosX`
- Modified parameters: `centromereSize`, `rulerIntervalMb`
- function `asymmetry` better dealing with unexpected cases
- function `armRatioCI` better dealing with unexpected cases

Shiny:

- New menu `stats` to export indices and marks’ stats
- Several parameters added
- Several examples added to presets

Bugs:

- Better dealing with horizontal karyotypes
- Fixed bug of displaying pos. of bands when `showMarkPos= TRUE`

param:

- `rulerIntervalMb`: Use data in millions
- `leftNotesUpPosX`: numeric, move up left notes to the right or left (x
  axis)
- `autoCenSize` boolean, when `TRUE` ignores `centromereSize`
- `centromereSize`: Apparent size of centromeres. Requires
  `autoCenSize = FALSE`
- `showBandList` in `posCalc` function to avoid adding mark names

# idiogramFISH 2.0.5

19-05-2021

Shiny:

- New dependence install options for the shiny app `runBoard()` function
- Better dealing with pandoc versions in shiny app

# idiogramFISH 2.0.4

23-04-2021

- Fixed bug of too many fonts (angles) when hundreds of marks for .svg
  for `labelOutwards=TRUE`
- jupyter notebooks for working locally or online in colab

# idiogramFISH 2.0.3

12-04-2021

- Downloading with `rentrez` package and plotting of chromosomes in
  shiny app
- Improvements in `genBankReadIF` function

param:

- `markN` numeric, vertices number for round corners of marks

# idiogramFISH 2.0.2

02-03-2021

- Added new style of mark `exProtein` (param: `startPos`, `pMarkFac`)
  and `inProtein` for circular plots also.
- Added new style of centromere `inProtein`
- better reading of genBank data by function `genBankReadIF`
- changes in dealing with `chrWidth` and other param. in circular plots
- Added parameters: `startPos`, `pMarkFac`,`cenFormat`,`cenFactor`
- Modified parameters: `bannedMarkNameAside`, `xModifier`, `roundedCen`

param:

- `bannedMarkNameAside` renamed to `bMarkNameAside`
- `startPos` numeric, factor to increase separation of `exProtein` marks
  to chromosome. Defaults to `0`
- `pMarkFac` numeric, fraction of chr. size for `exProtein` style marks.
  Defaults to `0.25`
- `xModifier` was modified, now quotient of `chrWidth`
- `roundedCen`: deprecated, see `cenFormat`
- `cenFormat`: character, when `"triangle"`, cen. has triangular aspect.
  When `"rounded"`, it has rounded aspect (Default). `"inProtein"` for
  using the mark with style of same name.
- `cenFactor`: numeric, modifies any cen. mark and cen. size. Defaults
  to `1`

Shiny App:

- Several changes in UI
- show/hide default values in code tab
- save data.frames as `.rds` too
- better update of log tab

# idiogramFISH 2.0.0

05-02-2021

- Added Shiny app, run with `runBoard()`
- Added option to add OTU name as `leftNotesUp` (`OTUasLeftNote`)

param:

- `anchorTextMoveParenX` numeric, for plots with `anchorTextMParental`
  move text in X axis. Defaults to `0`
- `anchorTextMoveParenY` numeric, for plots with `anchorTextMParental`
  move text in Y axis. Defaults to `0`
- `OTUasLeftNote` boolean, when `TRUE` adds OTU (species) name to the
  left-up

# idiogramFISH 1.16.8

30-11-2020

- Added rounded vertices for ggplot related functions
- Better parsing of F<sub>L</sub> (chr. name formatting) and similar
  text in notes
- Sorting of chr. by `chrNameUp` column added.

param:

- `leftNotesUpPosY` numeric, move up-left-notes `leftNotesUp` down or up
  (y axis)

# idiogramFISH 1.16.7

15-10-2020

- Helper functions for simple plots with ggplot `mapGGChr` and
  `mapGGChrMark`
- Better customization of `anchor`

param:

- `addMissingOTUBefore`: character, when you want to add space (ghost
  OTUs) before one or several OTUs, pass the names of OTUs after the
  desired space in a character vector
  i.e. `c("species one","species five")`
- `karAnchorRight`: character, OTUs’ add anchor to the right of this OTU
  (names of karyotypes). For `verticalPlot=FALSE`
- `anchorText`: character, text to add to `anchor` structure near
  symbol. See `anchor`. Defaults to `""`
- `anchorTextMParental`: character, designed to fill with a character
  object the space left of a missing parental in the `anchor` structure.
- `anchorTextMoveX`: numeric, for vertical plots with `anchorText` move
  text in X axis. Defaults to `0.5`
- `anchorTextMoveY`: numeric, for horizontal plots with `anchorText`
  move text in Y axis. Defaults to `1`
- `anchorLineLty`: numeric, type of line in `anchor`, corresponds to
  `lty`. Defaults to `1`

# idiogramFISH 1.16.6

16-09-2020

- Support for F<sup>+</sup> chr. types in citrus functions
- Horizontal arrange introduced `verticalPlot=FALSE`
- Horizontal anchor for karyotypes in `karAnchorLeft`
- Move all karyotypes `moveAllKarValueY`, `moveAllKarValueHor`
- Introducing more notes, over kar., `leftNotesUp`
- `dotsAsOval` convert dots style to one oval. not for circ. plots.
- `useOneDot` now for regular plots (non-circular) also.
- better plotting of white chr. in circular plots
- Fixed bug of dots size in circular plots
- Fixed bug of ruler of long arm.
- Fixed bug when chr. name F<sub>L</sub>
- Fixed bug when calculating mark size when NA - GISH
- Better plotting of white chr.
- square marks name splitting

param:

- `verticalPlot`: boolean, when `TRUE` karyotypes are plotted
  vertically, otherwise, horizontally. Defaults to `TRUE`
- `karSpaceHor`: numeric, separation among horizontal karyotypes. When
  `verticalPlot=FALSE`. Defaults to `0`
- `karAnchorLeft`: character, OTUs’ names of karyotypes to the right of
  your desired anchor. For `verticalPlot=FALSE`
- `moveAllKarValueHor`: numeric, similar to `mkhValue`, but affects all
  karyotypes.
- `moveAllKarValueY`: numeric, similar to `moveAllKarValueHor`, but
  affects y axis.
- `leftNotesUp`: data.frame, (to the left), similar to `leftNotes`, but
  intended for placement over kar.
- `leftNotesPosX`: (`0.5`) numeric, moves left notes in the x axis
- `notesPosX`: (`0.5`) numeric, moves right notes in the x axis
- `noteFont`: numeric `1` for normal, `2` for bold, `3` for italics, `4`
  for bold-italics. See `notes`
- `leftNoteFont`: numeric `1` for normal, `2` for bold, `3` for italics,
  `4` for bold-italics. See `leftNotes`
- `leftNoteFontUp`: numeric `1` for normal, `2` for bold, `3` for
  italics, `4` for bold-italics. See `leftNotesUp`
- `parseTypes`: boolean, parse in `notes` the Citrus chr. types names.
  Creates subindex pos. for FL.
- `parseStr2lang`: bolean, parse string in `notes` with function
  `str2lang(paste0("paste(",note,")") )` for ex:
  `"italic('C. sinensis'), ' Author'"`. See `notes`,
  `leftNotes`,`leftNotesUp`.
- `gishCenBorder`: boolean, when `TRUE`, cen. mark border color is the
  same as mark color, ignoring `colorBorderMark`. No default.
- `hideCenLines`: numeric, factor to multiply line width (lwd) used for
  covering cen. border, when `chrColor` is `white` or when
  `gishCenBorder=TRUE`
- `markNewLine`, character, character to split mark Names into different
  text lines. Applies to `square` marks. Defaults to `NA`
- `mylheight`, numeric, for `markNewLine!=NA`; is equivalent to
  `lheight` of `par`: “The line height multiplier. The height of a line
  of text (used to vertically space multi-line text) is found by
  multiplying the character height both by the current character
  expansion and by the line height multiplier.” Defaults to `0.7`.
- `bannedMarkNameAside`: boolean, when `TRUE` and `legend="inline"`,
  shows marks in `bannedMarkName` as `legend="aside"` would do. See
  `bannedMarkName`
- `forbiddenMark`: character, character string or vector with mark names
  to be removed from plot. Not the marks but the labels.
- `lwd.marks`: thickness of most marks. Except `cM` marks and centr.
  related marks. See `lwd.chr`, `lwd.cM`
- `dotsAsOval`: boolean, use oval instead of two dots in style of marks
  `dots`. Defaults to `FALSE`. See `useOneDot`. Not useful for
  `chromatids=TRUE` or `circularPlot=TRUE`

# idiogramFISH 1.16.1

29-07-2020

- `squareLeft` new style of mark. as `square` but with legend to the
  left when inline.
- minor ticks possible
- add mark % of chr. and position to plot.
- additional column `chrNameUp` for name over kar.
- show chr. sizes in μm and Mbp under kar.
- when `OTUfont=3` (italics), var. name present inside `'` is not shown
  in italics
- added anchor structure for progenies, see GISH
- `notesLeft` deprecated pass data.frame to `leftNotes`
- Better separation of groups
- F<sub>L</sub><sup>+</sup> chr. name now corrected in plot, previously
  F<sub>L</sub><sup>NA</sup>
- `ylabline` renamed to `xPosRulerTitle`

param:

- `groupSepar`: numeric, factor for affecting chr. spacing `chrSpacing`
  among groups. Defaults to `0.5`
- `useMinorTicks`: boolean, display minor ticks between labeled ticks in
  ruler. See `miniTickFactor`. Defaults to `FALSE`. (ticks without
  label)
- `miniTickFactor`: numeric, number of minor ticks for each labeled
  tick. See `useMinorTicks`. Defaults to `10`
- `xPosRulerTitle`: (`2.6`) Modifies the horizontal position of the
  title of rulers (Mb, etc). Moves to left from 1st chr. in `chrSpacing`
  times
- `yPosRulerTitle`: numeric, affects vertical position of ruler title.
  Defaults to `0`
- `markPer`: character, name of mark to calculate % of mark in chr. and
  add it to plot. See `perAsFraction`
- `perAsFraction`: boolean, when `TRUE` % is shown as fraction. Defaults
  to `FALSE`. See `markPer`
- `showMarkPos`: boolean, adds position of marks under karyotype
  (fraction 0-1) when `TRUE`. Defaults to `FALSE`
- `bToRemove`: character, bands to remove from calc. of pos.
- `chrSize` show chr. size under karyo.
- `chrNameUp` use col. of the same name to add secondary name over kar.
- `classMbName` “chromosome” name when in Mbp
- `classcMName` “chromosome” name when in cM
- `classChrName` “chromosome” name when in μm
- `classChrNameUp` “chromosome” name `chrNameUp`
- `classGroupName` name of title of groups
- `nsmall` digits for rounding of `chrSize`
- `chrSizeMbp` show chr. size Mbp requires col. `Mbp`
- `groupName`, hide or show group name
- `leftNotes`, similar to notes
- `leftNotesPos`, x
- `leftNotesPosY`, y
- `moveKarHor`, move kar. to right
- `mkhValue`, amount to move to right
- `anchor`, display anchor for moveKarHor OTUs
- `anchorVsizeF` factor to modify vertical size of anchor
- `moveAnchorV`, move anchor vertical portion
- `moveAnchorH`, move anchor horizontal portion
- `pchAnchor`, symbol in anchor
- `rulerPosMod` deprecated

# idiogramFISH 1.15.3

01-07-2020

main changes:

- Better plotting of GISH with chromatids
- change in messages when missing data
- chr. in groups are closer
- parsing of citrus names of chromosomes F<sub>L</sub><sup>+</sup> and
  F<sub>L</sub><sup>0</sup> automatic
- helper functions for plotting Citrus karyotypes:
  - `citrusSize`, `citrusMarkPos`, `markOverCMA`

param: (`plotIdiograms`)

- `efZero` threshold for checking if != 0
- `orderChr`, order of chr. Replaces `orderBySize` - deprecated. Values
  = `size`, `original`, `name`, `group`
- `orderBySize` - deprecated
- `notesLeft` note position to the left when `TRUE`
- `notesPosY` y axis modify notes position
- `chrIdPatternRem` regex pattern to remove from chr. names

# idiogramFISH 1.15.1

02-06-2020

- introducing ‘chromatids’
- new rounded style of centromere added (default).
- better naming of `w` position marks when `inline`
- changed logic of `cenStyle` coloring

param:

- `chromatids` show separated chromatids
- `holocenNotAsChromatids` do not use chromatids in holocen.
- `arrowsBothChrt` prints arrow marks in both chromatids
- `excHoloFrArrToSide` excludes holocentrics from arrowsToSide config.
- `xModifier` separation among chromatids
- `xModMonoHoloRate` shrink holocen. separ among chromatids with this
  quotient.
- `remSimiMarkLeg` remove “duplicated” name of labels when presence of
  pseudoduplicates arising from `pattern`
- `bannedMarkName` remove this mark name from labels (legends)
- `defCenStyleCol` color for external part of marks with `cenStyle`
- `roundedCen` rounded centromere
- `lwd.mimicCen` line width for `cenStyle` marks
- `squareness` new name for `roundness` (deprecated)

# idiogramFISH 1.14.11

23-04-2020

- `genBankReadIF` function, now allows duplicated field names
- `cMLeft` style of mark added
- `cM` and `cMLeft`styles are used as `inline` type of `legend` for
  arrows (`upArrow`,`downArrow`)
- A new column `protruding` can be added to `dfMarkColor` data.frame to
  define aspect of `cM` marks
- `namesToColumns` new function to avoid overlap of mark names, for
  holoc. and monoc.

params (`namesToColumns`):

- `marksDf` data.frame of marks
- `dfChrSize` data.frame, size of chr. same as plot.
- `markType` of type “downArrow”,“upArrow”,“cM”,“cMLeft”
- `amountofSpaces` numeric, number of spaces for each column
- `colNumber` numeric, number of columns
- `protruding` numeric, same as plot, equivalent to cM protruding
- `protrudingInt` numeric, spacing of columns in terms of width of chr.
  percent 1 = 100%.
- `circularPlot` same as plot
- `rotation` same as plot
- `defaultStyleMark` if some data in column style missing fill with this
  one
- `orderBySize` same as in plot.
- `halfModUp` when plotting several chromosomes in a circular plot,
  corrects for alignment problems of “upArrows”, “cM” labels.
- `halfModDown` when plotting several chromosomes in a circular plot,
  corrects for alignment problems of “downArrows”, “cMLeft” labels.
- `rotatMod` for circ. plots, when rotation diff. from `0`, corrects
  alignment of labels.

params:

- `cMBeginCenter` modifies start position of `cM` and `cMLeft` marks
- `arrowsToSide` arrows are plotted near chr. margin

# idiogramFISH 1.14.7

27-03-2020

- Compatibility with `rentrez` downloaded data
- Better reading of `join` from genBank data
- new styles of mark: `cenStyle` to add constrictions anywhere;
  `upArrow` (clockwise in circular plot); `downArrow` (anti-clockwise in
  circular plot)
- fixed bug when `legend="inline"` in circular plots

params:

- `rulerTitleSize`: Font size of units (title)
- `arrowhead`: proportion of head of arrow - length
- `shrinkArrow`: proportion to shrink body of arrow - width
- `arrowheadWidthShrink`: proportion to shrink arrowhead - width

# idiogramFISH 1.14.2

26-02-2020

- Introducing circular plots `circularPlot=TRUE` and other params. for
  circular plot
- function `genBankReadIF` to read plasmid or prokaryote data. Uses
  tidyr.  
- function `swapChrRegionDfSizeAndMarks` to swap arm size and marks
- tolerance when column `markSize` absent

params:

- `legendYcoord`: modify mark legend Y pos (for common plot also)
- `callPlot`: call plot.new or use your device (when FALSE)

params: (circularPlot=TRUE)

- `shrinkFactor`: size of chr. in fraction of circle
- `separFactor`: separ among kar.
- `labelSpacing`: among label and chr.
- `chrLabelSpacing`: chr. label space
- `OTUlabelSpacing`: OTU name space
- `radius`: radius
- `OTUsrt`: angle of OTU name text
- `OTUplacing`: add number and legend instead of OTU name
- `useOneDot`: one dot instead of two
- `circleCenter`: X coordinate
- `circleCenterY`: Y coordinate
- `OTULabelSpacerx`: modify OTU name pos.
- `OTULabelSpacery`: modify OTU name pos.
- `OTUcentered`: OTU name centered
- `OTUjustif`: OTU name justif.
- `OTUlegendHeight`: separ. of OTU names when `OTUplacing`

# idiogramFISH 1.13.8

05-02-2020

- Fixed bug when plotting several OTU with groups
- cen. marks allowed also when `centromereSize = 0`
- improvement in automatic scale of ruler.
- Added the “cM” style of mark, with custom `protruding`
- `centromereSize` is automatic (when absent), as well as
  `rulerInterval`

params:

- `lwd.cM`: thickness of cM marks
- `OTUfont`: style of font of OTU name
- `OTUfamily`: font family for OTU names
- `lwd.chr`: affects ruler too.
- `defaultFontFamily`: modify font of texts.
- Custom default style of mark with `defaultStyleMark`
- `fixCenBorder` affects cen. marks also.
- `chrBorderColor` for adding optionally chr. border color.
- `cenColor` defaults to `chrColor` now.
- `colorBorderMark` forces custom color in border of marks.
- `borderOfWhiteMarks`, if `TRUE`, when mark is white, its border is
  black.
- `ceilingFactor` number of significative digits to consider when
  rounding ruler max. value.
- `MbThreshold` created (substitutes `MbThresholds`)
- added option to modify ruler intervals for Mb, and cM independently
  with params: `rulerIntervalMb`, `rulerIntervalcM`
- other added parameters: `defaultStyleMark`, `protruding`,
  `ceilingFactor`, `rulerInterval`, `threshold`, `MbUnit`,
  `specialChrWidth`, `specialChrSpacing`, `specialOTUNames`,
  `specialyTitle`
- OTUs passed to `specialOTUNames`, can have special: `specialChrWidth`,
  `specialChrSpacing`, and `specialyTitle`. Useful for e.g. cM.
- Allowed customization of ruler (`ceilingFactor`, `rulerInterval`)
- Allowed custom ruler title `MbUnit`, `specialyTitle`, `yTitle`.
  `yTitle` is the common (micrometers). `specialyTitle` is for OTUs in
  `specialOTUNames` (e.g. “cM”), and `MbUnit` when data in millions and
  OTU is not in `specialOTUNames`

# idiogramFISH 1.12.1

06-01-2020

- Fixed bug of absence of chr. indices when monocen. and holocen.
  together
- Added functionality to print each index separately
- Added functionality to print groups below chr. name
- DOI added
- minor vignette corrections

# idiogramFISH 1.11.1

12 12 2019

- Added functionality for fixing y x aspect ratio (roundness proportion)
  using `asp = 1` only
- Use of `dotRoundCorr` discouraged, requires `useXYfactor = TRUE`
- Fixed misplacement of marks when `origin="t"` or
  `markDistType = "cen"`
- Added functionality for plotting karyotypes in micrometers and bases
  together, see monocen. vignette

# idiogramFISH 1.9.1

29 11 2019

- Fixed bug when centromere=0 when several karyotypes
- Added rounded vertices for `centromere > 0`
- Added functionality for plotting GISH.

# idiogramFISH 1.8.3

14 11 2019

- Fixed dependencies
- Fixed size of dots of legend

# idiogramFISH 1.8.1

29 10 2019

- Added parameters for adding notes to the right of karyotype.
- Improvement in messages when plotting.

# idiogramFISH 1.7.1

20 10 2019

- Cen. marks don’t need another data.frame. Can be present in main marks
  data.frame
- Allowed dup. names for not ordered chr. names (and no marks)

# idiogramFISH 1.6.3

13 10 2019

- More tolerance when allowing duplicated chr. names when no marks.
- Documentation changes, new examples.

# idiogramFISH 1.6.1

02 10 2019

- Added support to plot monocen. and holocen. together  
- Function `plotIdiogramsHolo` deprecated

# idiogramFISH 1.5.1

27 09 2019

- Added support to plot alongside phylogenies
- Allow some karyotypes to appear without indexes when error in long /
  short classif.
- Fixed bug in naming of OTUs.
- Vignettes corrections
- Fix references of packages in vignettes when package not installed.  
- Added support for vignettes in devel. in R-32 bits

# idiogramFISH 1.2.1

17 09 2019

- Fixed bug in armRatioCI that impacts all other functions.
- Added support for groups
- Added human karyotype
- Added rounded vertices when `centromereSize =0`
- You don’t have to use dfMarkColor data.frame, is not mandatory now.
- You can use (optionally) a character vector to pass colors.
- Package has default colors now.
