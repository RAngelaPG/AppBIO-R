observeEvent(input$confirmSvg, {
  removeModal()
  updateRadioButtons(session, inputId = "pngorsvg", selected = "svg")
  values[["pngorsvg"]] <- "svg"
  values[["decision"]] <- ". You did something not recommended, Now crashing?"
  values[["stop"]] <- FALSE
})

observeEvent(input$confirmPng, {
  removeModal()
  updateRadioButtons(session, inputId = "pngorsvg", selected = "png")
  values[["stop"]] <- FALSE
  values[["pngorsvg"]] <- "png"
  values[["decision"]] <- ""
})

observeEvent(input$pngorsvg,
  ignoreInit = TRUE,
  {
    if (input$pngorsvg == "svg" & values[["number"]] == 10) {
      showModal(modalDialog(
        title = "WARNING, Rstudio or browser crash risk",
        tagList(
          "You selected to display a large chr. in .svg format",
          br(),
          "you should leave .png in", tags$strong("Display"), "for thousands of marks",
          br(),
          "and", tags$strong("download plot"), "as .svg if desired"
        ),
        easyClose = TRUE,
        footer = list(
          actionButton("confirmSvg", "Choose .svg (not recommended)"),
          actionButton("confirmPng", "Leave .png as suggested")
        )
      ))
      values[["stop"]] <- TRUE
      values[["decision"]] <- ""
    } else {
      values[["pngorsvg"]] <- input$pngorsvg
      values[["stop"]] <- FALSE
      values[["decision"]] <- ""
    }
  }
)

filenameR <- eventReactive(
  {
    input$examIncrease
    input$examDecrease
    input$exampleButton
    values[["df1"]]
    values[["df1MStyle"]]
    values[["df1Mark"]]

    values[["notes"]]
    values[["leftNotes"]]
    values[["leftNotesUp"]]

    input$leftNoteFontUp
    input$leftNotesPosX
    input$leftNotesPosY
    input$leftNotesUpPosY

    input$leftNotesUpPosX
    input$notesPosY
    input$leftNoteFont
    input$noteFont
    input$leftNotesUpTextSize
    input$leftNotesTextSize

    input$verticalPlot
    input$karSpaceHor
    input$karAnchorLeft

    input$parseStr2lang
    input$moveAllKarValueHor
    input$moveAllKarValueY

    input$karHeight

    input$colorBorderMark
    input$alpha_val
    input$lwd.marks
    input$gishCenBorder
    input$hideCenLines

    input$karHeiSpace
    input$amoSepar
    input$karSepar
    input$cenFormat
    input$cenFactor

    input$autoCenSize
    input$centromereSize

    as.character(input$legend)
    values[["mycolors"]]
    values[["markPer"]]
    values[["bToRemove"]]
    values[["bannedMarkName"]]
    values[["specialOTUNames"]]
    values[["missOTUspacings"]]
    values[["addMissingOTUAfter"]]

    input$specialyTitle
    input$specialChrWidth
    input$specialChrSpacing

    input$origin
    input$OTUfamily
    input$classChrName
    input$classChrNameUp

    input$chrIdPatternRem
    input$xModMonoHoloRate

    input$yTitle
    input$nameChrIndexPos

    input$perAsFraction

    input$chrColor
    input$cenColor
    input$chrId
    input$chrWidth
    input$chrSpacing
    input$groupSepar
    input$squareness
    input$markDistType

    input$markN
    input$n
    input$orderChr
    input$useOneDot
    input$pMarkFac
    input$cMBeginCenter

    input$protruding
    input$arrowhead
    input$chrLabelSpacing
    input$labelSpacing

    input$rotation
    input$shrinkFactor
    input$separFactor
    input$radius
    input$circleCenter
    input$OTUsrt
    input$OTUjustif
    input$OTULabelSpacerx
    input$OTUlegendHeight
    input$OTUplacing
    input$morpho
    input$chrIndex
    input$karIndex
    input$karIndexPos
    input$chrSize
    input$showMarkPos
    input$indexIdTextSize
    input$chrSizeMbp
    input$miniTickFactor
    input$nsmall

    input$forbiddenMark

    input$chrNameUp
    input$useMinorTicks
    input$chromatids

    input$holocenNotAsChromatids

    input$xModifier
    input$circularPlot
    input$addOTUName
    as.numeric(input$OTUfont)
    as.character(input$moveKarHor)
    input$mkhValue
    input$anchor
    input$moveAnchorV
    input$moveAnchorH
    input$anchorVsizeF
    input$anchorHsizeF
    input$OTUasNote
    input$labelOutwards
    input$notesTextSize
    input$notesPosX
    input$ruler
    input$rulerPos
    input$threshold
    input$ceilingFactor

    input$rulerInterval
    input$rulerIntervalMb
    input$rulerIntervalcM

    input$ruler.tck
    input$collapseCen

    input$rulerNumberSize
    input$rulerNumberPos
    input$rulerTitleSize
    input$xPosRulerTitle
    input$legendWidth
    input$legendHeight
    input$pattern
    input$markNewLine
    input$classGroupName
    input$remSimiMarkLeg
    input$markLabelSize
    input$markLabelSpacer
    input$legendYcoord
    input$fixCenBorder
    input$bMarkNameAside
    input$chrBorderColor
    input$lwd.chr
    input$lwd.mimicCen
    input$lwd.cM
    input$distTextChr
    input$OTUTextSize
    input$xlimLeftMod
    input$xlimRightMod
    input$ylimBotMod
    input$ylimTopMod

    values[["pngorsvg"]]
  },
  {
    filename <- tempfile(fileext = ".txt")
    filenamePath <- normalizePath(filename, mustWork = F)
    filenamePath
  }
)

observeEvent(input$left, {
  updateNumericInput(session, "hwModifier", value = isolate(input$hwModifier) - 0.2)
})

observeEvent(input$right, {
  updateNumericInput(session, "hwModifier", value = isolate(input$hwModifier) + 0.2)
})

observeEvent(input$fontDecrease, {
  updateNumericInput(session, "rulerNumberSize", value = isolate(input$rulerNumberSize) - 0.2)

  updateNumericInput(session, "notesTextSize", value = isolate(input$notesTextSize) - 0.2)
  updateNumericInput(session, "leftNotesTextSize", value = isolate(input$leftNotesTextSize) - 0.2)
  updateNumericInput(session, "leftNotesUpTextSize", value = isolate(input$leftNotesUpTextSize) - 0.2)

  updateNumericInput(session, "OTUTextSize", value = isolate(input$OTUTextSize) - 0.2)
  updateNumericInput(session, "rulerTitleSize", value = isolate(input$rulerTitleSize) - 0.2)
  updateNumericInput(session, "markLabelSize", value = isolate(input$markLabelSize) - 0.2)
  updateNumericInput(session, "indexIdTextSize", value = isolate(input$indexIdTextSize) - 0.2)
})

observeEvent(input$fontIncrease, {
  updateNumericInput(session, "rulerNumberSize", value = isolate(input$rulerNumberSize) + 0.2)

  updateNumericInput(session, "notesTextSize", value = isolate(input$notesTextSize) + 0.2)
  updateNumericInput(session, "leftNotesTextSize", value = isolate(input$leftNotesTextSize) - 0.2)
  updateNumericInput(session, "leftNotesUpTextSize", value = isolate(input$leftNotesUpTextSize) - 0.2)

  updateNumericInput(session, "OTUTextSize", value = isolate(input$OTUTextSize) + 0.2)
  updateNumericInput(session, "rulerTitleSize", value = isolate(input$rulerTitleSize) + 0.2)
  updateNumericInput(session, "markLabelSize", value = isolate(input$markLabelSize) + 0.2)
  updateNumericInput(session, "indexIdTextSize", value = isolate(input$indexIdTextSize) + 0.2)
})

add0 <- "{"

addLibraryIdio <- "library(idiogramFISH)"

observeEvent(input$pngorsvgDown, {
  if (input$pngorsvgDown == "svg") {
    values[["imageType"]] <- "image/svg+xml"
  } else {
    values[["imageType"]] <- "image/png"
  }
})


observe({
  values[["mycolors"]] <- tryCatch(unlist(strsplit(input$mycolors, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["markPer"]] <- tryCatch(unlist(strsplit(input$markPer, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["bToRemove"]] <- tryCatch(unlist(strsplit(input$bToRemove, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["bannedMarkName"]] <- tryCatch(unlist(strsplit(input$bannedMarkName, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["specialOTUNames"]] <- tryCatch(unlist(strsplit(input$specialOTUNames, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["missOTUspacings"]] <- tryCatch(as.numeric(unlist(strsplit(input$missOTUspacings, ","))),
    error = function(e) {
      NULL
    }
  )

  values[["addMissingOTUAfter"]] <- tryCatch(unlist(strsplit(input$addMissingOTUAfter, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["mycolors2"]] <- tryCatch(unlist(strsplit(input$mycolors2, ",")),
    error = function(e) {
      NULL
    }
  )


  if (inherits(values[["df1"]], "data.frame") & !invalid(values[["df1"]])) {
    dfChrSizeName <- ifelse(input$saveType == "rds",
      input$chrfilename2,
      shQuote(paste0(isolate(input$chrfilename2), ".csv"))
    )

    chrRds <- paste0(input$chrfilename2, "  <- readRDS('", input$chrfilename2, ".rds')\n")
  } else {
    dfChrSizeName <- shQuote("")
    values[["df1"]] <- NULL
    chrRds <- ""
  }

  if (inherits(values[["notes"]], "data.frame") & !invalid(values[["notes"]])) {
    notesName <- tryCatch(ifelse(input$saveType == "rds",
      input$notesfilename2,
      shQuote(paste0(isolate(input$notesfilename2), ".csv"))
    ), error = function(e) {
      NULL
    })

    notesRds <- paste0(input$notesfilename2, " <- readRDS('", input$notesfilename2, ".rds')\n")
  } else {
    notesName <- shQuote("")
    values[["notes"]] <- NULL
    notesRds <- ""
  }

  if (inherits(values[["leftNotes"]], "data.frame") & !invalid(values[["leftNotes"]])) {
    leftNotesName <- tryCatch(ifelse(input$saveType == "rds",
      input$leftNotesfilename2,
      shQuote(paste0(isolate(input$leftNotesfilename2), ".csv"))
    ), error = function(e) {
      NULL
    })
    leftNotesRds <- paste0(input$leftNotesfilename2, " <- readRDS('", input$leftNotesfilename2, ".rds')\n")
  } else {
    leftNotesName <- shQuote("")
    values[["leftNotes"]] <- NULL
    leftNotesRds <- ""
  }

  if (inherits(values[["leftNotesUp"]], "data.frame") & !invalid(values[["leftNotesUp"]])) {
    leftNotesUpName <- tryCatch(ifelse(input$saveType == "rds",
      input$leftNotesUpfilename2,
      shQuote(paste0(isolate(input$leftNotesUpfilename2), ".csv"))
    ), error = function(e) {
      NULL
    })

    leftNotesUpRds <- paste0(input$leftNotesUpfilename2, " <- readRDS('", input$leftNotesUpfilename2, ".rds')\n")
  } else {
    leftNotesUpName <- shQuote("")
    values[["leftNotesUp"]] <- NULL
    leftNotesUpRds <- ""
  }

  if (inherits(values[["df1Mark"]], "data.frame") & !invalid(values[["df1Mark"]])) {
    dfMarkPosName <- tryCatch(ifelse(input$saveType == "rds",
      input$markfilename2,
      shQuote(paste0(isolate(input$markfilename2), ".csv"))
    ), error = function(e) {
      NULL
    })

    markRds <- paste0(input$markfilename2, " <- readRDS('", input$markfilename2, ".rds')\n")
  } else {
    dfMarkPosName <- shQuote("")
    values[["df1Mark"]] <- NULL
    markRds <- ""
  }

  if (inherits(values[["df1MStyle"]], "data.frame") & !invalid(values[["df1MStyle"]])) {
    dfMarkColorName <- tryCatch(ifelse(input$saveType == "rds",
      input$MStylefilename2,
      shQuote(paste0(isolate(input$MStylefilename2), ".csv"))
    ), error = function(e) {
      NULL
    })

    mstyleRds <- paste0(input$MStylefilename2, " <- readRDS('", input$MStylefilename2, ".rds')\n")
  } else {
    dfMarkColorName <- shQuote("")
    values[["df1MStyle"]] <- NULL
    mstyleRds <- ""
  }

  mc <- match.call(
    plotIdiograms,
    call("plotIdiograms",
      dfChrSize = dfChrSizeName,
      notes = notesName,
      leftNotes = leftNotesName,
      leftNotesUp = leftNotesUpName,
      dfMarkPos = dfMarkPosName,
      dfMarkColor = dfMarkColorName,
      colorBorderMark = shQuote(input$colorBorderMark),
      alpha_val = input$alpha_val,
      lwd.marks = input$lwd.marks,
      gishCenBorder = input$gishCenBorder,
      hideCenLines = input$hideCenLines,
      leftNoteFontUp = input$leftNoteFontUp,
      leftNotesPosX = input$leftNotesPosX,
      leftNotesPosY = input$leftNotesPosY,
      leftNotesUpPosY = input$leftNotesUpPosY,
      leftNotesUpPosX = input$leftNotesUpPosX,
      notesPosY = input$notesPosY,
      leftNoteFont = input$leftNoteFont,
      noteFont = input$noteFont,
      leftNotesUpTextSize = input$leftNotesUpTextSize,
      leftNotesTextSize = input$leftNotesTextSize,
      verticalPlot = input$verticalPlot,
      karSpaceHor = input$karSpaceHor,
      karAnchorLeft = shQuote(input$karAnchorLeft),
      parseStr2lang = input$parseStr2lang,
      moveAllKarValueHor = input$moveAllKarValueHor,
      moveAllKarValueY = input$moveAllKarValueY,
      karHeight = input$karHeight,
      karHeiSpace = input$karHeiSpace,
      amoSepar = input$amoSepar,
      karSepar = input$karSepar,
      legend = shQuote(input$legend),
      mycolors = if (length(values[["mycolors"]]) == 0) {
        shQuote("")
      } else if (length(values[["mycolors"]]) == 1) {
        shQuote(values[["mycolors"]])
      } else {
        values[["mycolors"]]
      },
      markPer = if (length(values[["markPer"]]) == 0) {
        shQuote("")
      } else if (length(values[["markPer"]]) == 1) {
        shQuote(values[["markPer"]])
      } else {
        values[["markPer"]]
      },
      bToRemove = if (length(values[["bToRemove"]]) == 0) {
        shQuote("")
      } else if (length(values[["bToRemove"]]) == 1) {
        shQuote(values[["bToRemove"]])
      } else {
        values[["bToRemove"]]
      },
      bannedMarkName = if (length(values[["bannedMarkName"]]) == 0) {
        shQuote("")
      } else if (length(values[["bannedMarkName"]]) == 1) {
        shQuote(values[["bannedMarkName"]])
      } else {
        values[["bannedMarkName"]]
      },
      specialOTUNames = if (length(values[["specialOTUNames"]]) == 0) {
        shQuote("")
      } else if (length(values[["specialOTUNames"]]) == 1) {
        shQuote(values[["specialOTUNames"]])
      } else {
        values[["specialOTUNames"]]
      },
      missOTUspacings = if (length(values[["missOTUspacings"]]) == 0) {
        ""
      } else {
        values[["missOTUspacings"]]
      },
      addMissingOTUAfter = if (length(values[["addMissingOTUAfter"]]) == 0) {
        shQuote("")
      } else if (length(values[["addMissingOTUAfter"]]) == 1) {
        shQuote(values[["addMissingOTUAfter"]])
      } else {
        values[["addMissingOTUAfter"]]
      },
      specialyTitle = shQuote(input$specialyTitle),
      specialChrWidth = input$specialChrWidth,
      specialChrSpacing = input$specialChrSpacing,
      origin = shQuote(input$origin),
      OTUfamily = shQuote(input$OTUfamily),
      classChrName = shQuote(input$classChrName),
      classChrNameUp = shQuote(input$classChrNameUp),
      chrIdPatternRem = shQuote(input$chrIdPatternRem),
      xModMonoHoloRate = input$xModMonoHoloRate,
      yTitle = shQuote(input$yTitle),
      nameChrIndexPos = input$nameChrIndexPos,
      perAsFraction = input$perAsFraction,
      cenFormat = shQuote(input$cenFormat),
      cenFactor = input$cenFactor,
      autoCenSize = input$autoCenSize,
      centromereSize = input$centromereSize,
      chrColor = shQuote(input$chrColor),
      cenColor = shQuote(input$cenColor),
      chrId = shQuote(input$chrId),
      chrWidth = input$chrWidth,
      chrSpacing = input$chrSpacing,
      groupSepar = input$groupSepar,
      squareness = input$squareness,
      markDistType = shQuote(input$markDistType),
      markN = input$markN,
      n = input$n,
      orderChr = shQuote(input$orderChr),
      useOneDot = input$useOneDot,
      pMarkFac = input$pMarkFac,
      protruding = input$protruding,
      arrowhead = input$arrowhead,
      cMBeginCenter = input$cMBeginCenter,
      chrLabelSpacing = input$chrLabelSpacing,
      labelSpacing = input$labelSpacing,
      rotation = input$rotation,
      shrinkFactor = input$shrinkFactor,
      separFactor = input$separFactor,
      radius = input$radius,
      circleCenter = input$circleCenter,
      OTUsrt = input$OTUsrt,
      OTUjustif = input$OTUjustif,
      OTULabelSpacerx = input$OTULabelSpacerx,
      OTUlegendHeight = input$OTUlegendHeight,
      OTUplacing = shQuote(input$OTUplacing),
      morpho = shQuote(input$morpho),
      chrIndex = shQuote(input$chrIndex),
      karIndex = input$karIndex,
      karIndexPos = input$karIndexPos,
      chrSize = input$chrSize,
      miniTickFactor = input$miniTickFactor,
      nsmall = input$nsmall,
      showMarkPos = input$showMarkPos,
      indexIdTextSize = input$indexIdTextSize,
      chrSizeMbp = input$chrSizeMbp,
      useMinorTicks = input$useMinorTicks,
      chrNameUp = input$chrNameUp,
      chromatids = input$chromatids,
      holocenNotAsChromatids = input$holocenNotAsChromatids,
      xModifier = input$xModifier,
      circularPlot = input$circularPlot,
      addOTUName = input$addOTUName,
      OTUfont = as.numeric(input$OTUfont),
      moveKarHor = shQuote(input$moveKarHor),
      mkhValue = input$mkhValue,
      anchor = input$anchor,
      moveAnchorV = input$moveAnchorV,
      moveAnchorH = input$moveAnchorH,
      anchorVsizeF = input$anchorVsizeF,
      anchorHsizeF = input$anchorHsizeF,
      OTUasNote = input$OTUasNote,
      labelOutwards = input$labelOutwards,
      notesTextSize = input$notesTextSize,
      notesPosX = input$notesPosX,
      ruler = input$ruler,
      rulerPos = input$rulerPos,
      threshold = input$threshold,
      ceilingFactor = input$ceilingFactor,
      rulerInterval = input$rulerInterval,
      rulerIntervalMb = input$rulerIntervalMb,
      rulerIntervalcM = input$rulerIntervalcM,
      ruler.tck = input$ruler.tck,
      collapseCen = input$collapseCen,
      rulerNumberSize = input$rulerNumberSize,
      rulerNumberPos = input$rulerNumberPos,
      rulerTitleSize = input$rulerTitleSize,
      xPosRulerTitle = input$xPosRulerTitle,
      legendWidth = input$legendWidth,
      legendHeight = input$legendHeight,
      pattern = shQuote(input$pattern),
      forbiddenMark = shQuote(input$forbiddenMark),
      markNewLine = shQuote(input$markNewLine),
      classGroupName = shQuote(input$classGroupName),
      remSimiMarkLeg = input$remSimiMarkLeg,
      markLabelSize = input$markLabelSize,
      markLabelSpacer = input$markLabelSpacer,
      legendYcoord = input$legendYcoord,
      fixCenBorder = input$fixCenBorder,
      bMarkNameAside = input$bMarkNameAside,
      chrBorderColor = shQuote(input$chrBorderColor),
      lwd.chr = input$lwd.chr,
      lwd.mimicCen = input$lwd.mimicCen,
      lwd.cM = input$lwd.cM,
      distTextChr = input$distTextChr,
      OTUTextSize = input$OTUTextSize,
      xlimLeftMod = input$xlimLeftMod,
      xlimRightMod = input$xlimRightMod,
      ylimBotMod = input$ylimBotMod,
      ylimTopMod = input$ylimTopMod
    )
  )

  mclist <- as.list(mc)
  mclist[1] <- NULL

  if (!invalid(!as.logical(input$keepDefault))) {
    if (!as.logical(input$keepDefault)) {
      param_only <- gsub("Default", "", names(paramValues))

      for (i in 1:length(param_only)) {
        myl <- unlist(mclist[which(names(mclist) == param_only[i])])

        myDef <- paramValues[i]

        if (!invalid(myl)) {
          if (!invalid(myDef[[1]])) {
            if (length(setequal(myl, myDef[[1]]))) {
              if (setequal(myl, myDef[[1]])) {
                mclist[which(names(mclist) == param_only[i])] <- NA
              }
            }
            if (length(setequal(myl, paste0("'", myDef[[1]], "'")))) {
              if (setequal(myl, paste0("'", myDef[[1]], "'"))) {
                mclist[which(names(mclist) == param_only[i])] <- NA
              }
            }
          }
          if (length(setequal(myl, ""))) {
            if (setequal(myl, "")) {
              mclist[which(names(mclist) == param_only[i])] <- NA
            }
          }
          if (length(setequal(myl, "''"))) {
            if (setequal(myl, "''")) {
              mclist[which(names(mclist) == param_only[i])] <- NA
            }
          }
        }
      }

      mclist <- mclist[which(sapply(mclist, unlist) != "''")]
      mclist <- mclist[which(!is.na(sapply(mclist, unlist)))]
    }
  }

  seq <- paste(names(mclist), mclist, sep = "=", collapse = ",\n")

  validate(
    need(try("dfChrSize" %in% names(mclist)), "")
  )

  if (!"dfChrSize" %in% names(mclist)) {
    seq <- character()
    mclist <- NULL
  } else if (!is.na(mclist$dfChrSize)) {
    if (mclist$dfChrSize == "''") {
      seq <- character()
      mclist <- NULL
    }
  } else {
    seq <- character()
    mclist <- NULL
  }

  if (!is.null(mclist)) {
    if ("dfChrSize" %in% names(mclist)) {
      if (inherits(values[["df1"]], "data.frame")) {
        mclist$dfChrSize <- values[["df1"]]
      }
    }
    if ("dfMarkPos" %in% names(mclist)) {
      if (inherits(values[["df1Mark"]], "data.frame")) {
        mclist$dfMarkPos <- values[["df1Mark"]]
      }
    }
    if ("dfMarkColor" %in% names(mclist)) {
      if (inherits(values[["df1MStyle"]], "data.frame")) {
        mclist$dfMarkColor <- values[["df1MStyle"]]
      }
    }
    if ("notes" %in% names(mclist)) {
      if (inherits(values[["notes"]], "data.frame")) {
        mclist$notes <- values[["notes"]]
      }
    }
    if ("leftNotes" %in% names(mclist)) {
      if (inherits(values[["leftNotes"]], "data.frame")) {
        mclist$leftNotes <- values[["leftNotes"]]
      }
    }
    if ("leftNotesUp" %in% names(mclist)) {
      if (inherits(values[["leftNotesUp"]], "data.frame")) {
        mclist$leftNotesUp <- values[["leftNotesUp"]]
      }
    }

    for (i in 1:length(mclist)) {
      if (is.character(mclist[[i]])) {
        mclist[[i]] <- sub("^'(.*)'$", "\\1", mclist[[i]])
      }
    }
  }

  addLibrarySvg <- "#install.packages('svglite')"
  addLibrarySvg <- ifelse(values[["pngorsvg"]] == "svg", addLibrarySvg, "")

  addLibrarySvg2 <- "library(svglite)"
  addLibrarySvg2 <- ifelse(values[["pngorsvg"]] == "svg", addLibrarySvg2, "")

  addSvgDev <- paste0(
    'svg("dfOfChrSize.svg",width=',
    values[["mysvgwidth"]],
    ", height=",
    values[["mysvgheight"]], ")"
  )

  addPngDev <- paste0(
    'png("dfOfChrSize.png",width=',
    values[["mysvgwidth"]] * 80,
    ", height=",
    values[["mysvgheight"]] * 80, ")"
  )

  addDev <- ifelse(values[["pngorsvg"]] == "svg", addSvgDev, addPngDev)

  mclist$pngorsvg <- input$pngorsvg
  mclist$heiFactor <- input$heiFactor
  mclist$widFactor <- input$widFactor

  rdsAdd <- ifelse(input$saveType == "rds", chrRds, "")
  notesAdd <- ifelse(input$saveType == "rds", notesRds, "")
  leftNotesAdd <- ifelse(input$saveType == "rds", leftNotesRds, "")
  leftNotesUpAdd <- ifelse(input$saveType == "rds", leftNotesUpRds, "")
  mrdsAdd <- ifelse(input$saveType == "rds", markRds, "")
  msrdsAdd <- ifelse(input$saveType == "rds", mstyleRds, "")

  block <- ifelse(input$asFile, "", "#")

  strFun <- paste0(
    add0, "\n",
    addLibrarySvg, "\n",
    block,
    addLibrarySvg2, "\n",
    addLibraryIdio, "\n",
    rdsAdd,
    mrdsAdd,
    msrdsAdd,
    "\n",
    notesAdd,
    leftNotesAdd,
    leftNotesUpAdd,
    block,
    ifelse(length(values[["mysvgwidth"]]), addDev, ""),
    "\n",
    ifelse(length(seq), paste0("plotIdiograms(\n", seq, "\n)"), ""),
    "\n",
    block,
    ifelse(length(values[["mysvgwidth"]]), "dev.off()", ""),
    "\n}"
  )

  strFun <- gsub("\n\n\n", "\n", strFun)
  strFun <- gsub("\n\n", "\n", strFun)

  spl1 <- unlist(strsplit(strFun, "\n"))
  spl2 <- character()

  for (s in spl1) {
    if (grepl("readRDS", s)) {
      to_s1 <- sub("(.*)<-.*", "\\1", s)
      to_s2 <- sub(".*(<-.*)", "\\1", s)
      to_s3 <- sprintf("%-*s%s%*s", 0, "", to_s1, ifelse((25 - nchar(to_s1)) < 1, 1, (25 - nchar(to_s1))), "")
      spl <- paste(to_s3, to_s2)
      spl2 <- c(spl2, spl)
    } else {
      spl2 <- c(spl2, s)
    }
  }

  strFun <- paste(spl2, collapse = "\n")

  if (length(input$keepDesc) > 0) {
    if (input$keepDesc) {
      spl1 <- unlist(strsplit(strFun, "\n"))
      spl2 <- character()

      for (s in spl1) {
        to_s <- sub("(.*)=.*", "\\1", s)
        to_ss <- sprintf("%-*s%s%*s", 0, "", s, ifelse((30 - nchar(s)) < 1, 1, (30 - nchar(s))), "")
        desc <- tryCatch(get(paste0(to_s, "Desc")), error = function(e) {
          ""
        })
        spl2 <- c(spl2, ifelse(desc == "", s, paste0(to_ss, "#", sub(paste0("`", to_s, "`:( \\([[:alnum:].`\u0022-]+\\))?"), "", desc))))
      }
      strFun <- paste(spl2, collapse = "\n")
    }
  }

  values[["strFun"]] <- strFun
  values[["presets"]] <- mclist
  # gc()
}) # observe 257
