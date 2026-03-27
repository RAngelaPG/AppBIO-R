observeEvent(input$nucPreset, {
  sel <- input$nucPreset
  updateRadioSubgroup(session, "exampleId", "sub", selected = sel, inline = TRUE)
})

observeEvent(input$exampleId, {
  sel <- input$exampleId
  sel2 <- as.numeric(gsub(".*-", "", sel))
  if (sel2 %in% 9:10) {
    updateRadioButtons(session, "nucPreset",
      selected = sel2, inline = TRUE
    )
  }
})

updateRadioSubgroup <- function(session, inputId, id = "sub", inline, selected, ...) {
  value <- paste0(id, "-", selected)
  updateRadioButtons(session, inputId, label = NULL, choices = NULL, inline = inline, selected = value)
}

radioGroupContainer <- function(inputId, ...) {
  class <- "form-group shiny-input-radiogroup shiny-input-container"
  div(id = inputId, class = class, ...)
}

radioSubgroup <- function(inputId, label, choiceNames, choiceValues, inline = TRUE, first = FALSE, id = "sub") {
  choiceValues <- paste0(id, "-", choiceValues)
  choices <- setNames(choiceValues, choiceNames)

  if (first == FALSE) {
    rb <- radioButtons(inputId, label, choices, selected = character(0), inline = inline)
  } else {
    rb <- radioButtons(inputId, label, choices, selected = choices[1], inline = inline)
  }
  rb$children
}

exampleId <- reactive({
  req(input$exampleId)
  parts <- unlist(strsplit(input$exampleId, "-"))
  list(
    id = parts[1],
    value = as.numeric(parts[2])
  )
})

#
#   param presets
#

observeEvent(
  {
    input$examIncrease
  },
  ignoreInit = TRUE,
  {
    if (values[["number"]] == 0 | values[["number"]] >= maxEx) {
      values[["number"]] <- 1
      values[["canButton"]] <- TRUE
    } else if (values[["number"]] <= maxEx - 1) {
      idx <- grep(paste0("\\<", exampleVec[exampleVec == values[["number"]]], "\\>"), exampleVec) + 1
      values[["number"]] <- exampleVec[idx]
      values[["canButton"]] <- TRUE
    }

    if (values[["number"]] %in% 9:10) {
      CurrentM$menu <- menulist[2]
      updateTabItems(session, "tabs", menulist[2])
      updateTabItems(session, "tabsetpanel4", tablist4[1])
    }
    values[["dfList"]] <- list()
    values[["button3ab"]] <- NA
    values[["entrez_summary1"]] <- NA
    values[["entrez_titles"]] <- NA
    values[["entrez_selected"]] <- NA
    values[["fetch_listAll"]] <- NA
    values[["authors"]] <- ""
    values[["entrez_search1"]] <- NA
    values[["fetch_list"]] <- NA
    values[["titles_number"]] <- NA
  }
)

observeEvent(
  {
    input$examDecrease
  },
  ignoreInit = TRUE,
  {
    if (values[["number"]] >= 2 & values[["number"]] <= maxEx) {
      idx <- grep(paste0("\\<", exampleVec[exampleVec == values[["number"]]], "\\>"), exampleVec) - 1
      values[["number"]] <- exampleVec[idx]
      values[["canButton"]] <- TRUE
    } else {
      values[["number"]] <- maxEx
      values[["canButton"]] <- TRUE
    }

    if (values[["number"]] %in% 9:10) {
      CurrentM$menu <- menulist[2]
      updateTabItems(session, "tabs", menulist[2])
      updateTabItems(session, "tabsetpanel4", tablist4[1])
    }
    values[["dfList"]] <- list()
    values[["button3ab"]] <- NA
    values[["entrez_summary1"]] <- NA
    values[["entrez_titles"]] <- NA
    values[["entrez_selected"]] <- NA
    values[["fetch_listAll"]] <- NA
    values[["authors"]] <- ""
    values[["entrez_search1"]] <- NA
    values[["fetch_list"]] <- NA
    values[["titles_number"]] <- NA
  }
)

observeEvent(
  {
    input$upPresetButton
  },
  ignoreInit = TRUE,
  {
    # values[["go"]]<-FALSE
    values[["canButton"]] <- FALSE
    updateCheckboxInput(session, "fileInsteadNcbi", value = FALSE)

    values[["number"]] <- input$upPreset
    values[["dfList"]] <- list()
    values[["button3ab"]] <- NA
    values[["entrez_summary1"]] <- NA
    values[["entrez_titles"]] <- NA
    values[["entrez_selected"]] <- NA
    values[["fetch_listAll"]] <- NA
    values[["authors"]] <- ""
    values[["entrez_search1"]] <- NA
    values[["fetch_list"]] <- NA
    values[["titles_number"]] <- NA

    values[["df1"]] <- values[["paramVec"]]$dfChrSizeVec[as.numeric(values[["number"]])][[1]]
    values[["df1Mark"]] <- values[["paramVec"]]$dfMarkPosVec[as.numeric(values[["number"]])][[1]]
    values[["df1MStyle"]] <- values[["paramVec"]]$dfMarkColorVec[as.numeric(values[["number"]])][[1]]
    values[["notes"]] <- values[["paramVec"]]$notesVec[as.numeric(values[["number"]])][[1]]
    values[["leftNotes"]] <- values[["paramVec"]]$leftNotesVec[as.numeric(values[["number"]])][[1]]
    values[["leftNotesUp"]] <- values[["paramVec"]]$leftNotesUpVec[as.numeric(values[["number"]])][[1]]
  }
)

observeEvent(
  {
    input$exampleButton
  },
  ignoreInit = TRUE,
  {
    values[["canButton"]] <- TRUE
    updateCheckboxInput(session, "fileInsteadNcbi", value = FALSE)

    values[["number"]] <- exampleId()$value # input$exampleId
    values[["dfList"]] <- list()
    values[["button3ab"]] <- NA
    values[["entrez_summary1"]] <- NA
    values[["entrez_titles"]] <- NA
    values[["entrez_selected"]] <- NA
    values[["fetch_listAll"]] <- NA
    values[["authors"]] <- ""
    values[["entrez_search1"]] <- NA
    values[["fetch_list"]] <- NA
    values[["titles_number"]] <- NA
  }
)

observeEvent(
  {
    input$nucPresetButton
  },
  ignoreInit = TRUE,
  {
    values[["canButton"]] <- TRUE
    updateCheckboxInput(session, "fileInsteadNcbi", value = FALSE)

    values[["number"]] <- input$nucPreset
    values[["dfList"]] <- list()

    values[["button3ab"]] <- NA
    values[["entrez_summary1"]] <- NA
    values[["entrez_titles"]] <- NA
    values[["entrez_selected"]] <- NA
    values[["fetch_listAll"]] <- NA
    values[["authors"]] <- ""
    values[["entrez_search1"]] <- NA
    values[["fetch_list"]] <- NA
    values[["titles_number"]] <- NA
  }
)

observeEvent(
  numberRe(),
  ignoreInit = TRUE,
  {
    validate(need(try(numberRe() > 0), "not ready"))

    sel <- numberRe()
    if (sel %in% 1:maxEx) {
      updateRadioSubgroup(session, "exampleId", "sub", selected = sel, inline = TRUE)
    }
    if (sel %in% 9:10) {
      updateRadioButtons(session, "nucPreset", selected = sel, inline = TRUE)
    }
    updateNumericInput(session, "karHeight", value = values[["paramVec"]]$karHeightVec[as.numeric(sel)])
    updateNumericInput(session, "karHeiSpace", value = values[["paramVec"]]$karHeiSpaceVec[as.numeric(sel)])
    updateNumericInput(session, "amoSepar", value = values[["paramVec"]]$amoSeparVec[as.numeric(sel)])
    updateNumericInput(session, "karSepar", value = values[["paramVec"]]$karSeparVec[as.numeric(sel)])
    updateNumericInput(session, "chrWidth", value = values[["paramVec"]]$chrWidthVec[as.numeric(sel)])
    updateRadioButtons(session, "chrId", selected = values[["paramVec"]]$chrIdVec[as.numeric(sel)])

    updateNumericInput(session, "squareness", value = values[["paramVec"]]$squarenessVec[as.numeric(sel)])
    updateNumericInput(session, "markN", value = values[["paramVec"]]$markNVec[as.numeric(sel)])
    updateNumericInput(session, "n", value = values[["paramVec"]]$nVec[as.numeric(sel)])

    updateRadioButtons(session, "orderChr", selected = values[["paramVec"]]$orderChrVec[as.numeric(sel)])
    updateRadioButtons(session, "markDistType", selected = values[["paramVec"]]$markDistTypeVec[as.numeric(sel)])

    updateCheckboxInput(session, "useOneDot", value = values[["paramVec"]]$useOneDotVec[as.numeric(sel)])
    updateTextInput(session, "bannedMarkName", value = values[["paramVec"]]$bannedMarkNameVec[as.numeric(sel)])
    updateTextInput(session, "specialOTUNames", value = values[["paramVec"]]$specialOTUNamesVec[as.numeric(sel)])
    updateTextInput(session, "addMissingOTUAfter", value = values[["paramVec"]]$addMissingOTUAfterVec[as.numeric(sel)])
    updateTextInput(session, "missOTUspacings", value = values[["paramVec"]]$missOTUspacingsVec[as.numeric(sel)])

    updateCheckboxInput(session, "origin", value = values[["paramVec"]]$originVec[as.numeric(sel)])
    updateTextInput(session, "OTUfamily", value = values[["paramVec"]]$OTUfamilyVec[as.numeric(sel)])

    updateTextInput(session, "classChrName", value = values[["paramVec"]]$classChrNameVec[as.numeric(sel)])
    updateTextInput(session, "classChrNameUp", value = values[["paramVec"]]$classChrNameUpVec[as.numeric(sel)])

    updateTextInput(session, "chrIdPatternRem", value = values[["paramVec"]]$chrIdPatternRemVec[as.numeric(sel)])

    updateNumericInput(session, "xModMonoHoloRate", value = values[["paramVec"]]$xModMonoHoloRateVec[as.numeric(sel)])

    updateTextInput(session, "specialyTitle", value = values[["paramVec"]]$specialyTitleVec[as.numeric(sel)])
    updateNumericInput(session, "specialChrWidth", value = values[["paramVec"]]$specialChrWidthVec[as.numeric(sel)])
    updateNumericInput(session, "specialChrSpacing", value = values[["paramVec"]]$specialChrSpacingVec[as.numeric(sel)])

    updateTextInput(session, "yTitle", value = values[["paramVec"]]$yTitleVec[as.numeric(sel)])
    updateNumericInput(session, "nameChrIndexPos", value = values[["paramVec"]]$nameChrIndexPosVec[as.numeric(sel)])

    updateCheckboxInput(session, "cMBeginCenter", value = values[["paramVec"]]$cMBeginCenterVec[as.numeric(sel)])
    updateNumericInput(session, "pMarkFac", value = values[["paramVec"]]$pMarkFacVec[as.numeric(sel)])

    updateNumericInput(session, "protruding", value = values[["paramVec"]]$protrudingVec[as.numeric(sel)])
    updateNumericInput(session, "arrowhead", value = values[["paramVec"]]$arrowheadVec[as.numeric(sel)])

    # circularPlot
    updateNumericInput(session, "chrLabelSpacing", value = values[["paramVec"]]$chrLabelSpacingVec[as.numeric(sel)])
    updateNumericInput(session, "labelSpacing", value = values[["paramVec"]]$labelSpacingVec[as.numeric(sel)])

    updateNumericInput(session, "rotation", value = values[["paramVec"]]$rotationVec[as.numeric(sel)])
    updateNumericInput(session, "shrinkFactor", value = values[["paramVec"]]$shrinkFactorVec[as.numeric(sel)])
    updateNumericInput(session, "separFactor", value = values[["paramVec"]]$separFactorVec[as.numeric(sel)])

    updateNumericInput(session, "radius", value = values[["paramVec"]]$radiusVec[as.numeric(sel)])
    updateNumericInput(session, "circleCenter", value = values[["paramVec"]]$circleCenterVec[as.numeric(sel)])

    updateNumericInput(session, "OTUsrt", value = values[["paramVec"]]$OTUsrtVec[as.numeric(sel)])
    updateNumericInput(session, "OTUjustif", value = values[["paramVec"]]$OTUjustifVec[as.numeric(sel)])
    updateNumericInput(session, "OTULabelSpacerx", value = values[["paramVec"]]$OTULabelSpacerxVec[as.numeric(sel)])
    updateNumericInput(session, "OTUlegendHeight", value = values[["paramVec"]]$OTUlegendHeightVec[as.numeric(sel)])
    updateRadioButtons(session, "OTUplacing", selected = values[["paramVec"]]$OTUplacingVec[as.numeric(sel)])

    # general
    updateNumericInput(session, "chrSpacing", value = as.numeric(values[["paramVec"]]$chrSpacingVec[as.numeric(sel)]))
    updateNumericInput(session, "groupSepar", value = as.numeric(values[["paramVec"]]$groupSeparVec[as.numeric(sel)]))

    updateRadioButtons(session, "morpho", selected = values[["paramVec"]]$morphoVec[as.numeric(sel)])
    updateTextInput(session, "chrColor", value = values[["paramVec"]]$chrColorVec[as.numeric(sel)])
    updateTextInput(session, "cenColor", value = values[["paramVec"]]$cenColorVec[as.numeric(sel)])

    updateRadioButtons(session, "chrIndex", selected = values[["paramVec"]]$chrIndexVec[as.numeric(sel)])
    updateCheckboxInput(session, "karIndex", value = values[["paramVec"]]$karIndexVec[as.numeric(sel)])
    updateNumericInput(session, "karIndexPos", value = values[["paramVec"]]$karIndexPosVec[as.numeric(sel)])

    updateCheckboxInput(session, "chrSize", value = as.logical(values[["paramVec"]]$chrSizeVec[as.numeric(sel)]))
    updateCheckboxInput(session, "showMarkPos", value = as.logical(values[["paramVec"]]$showMarkPosVec[as.numeric(sel)]))
    updateCheckboxInput(session, "useMinorTicks", value = as.logical(values[["paramVec"]]$useMinorTicksVec[as.numeric(sel)]))

    updateNumericInput(session, "indexIdTextSize", value = as.numeric((values[["paramVec"]]$indexIdTextSizeVec[as.numeric(sel)])))
    updateNumericInput(session, "miniTickFactor", value = as.numeric(values[["paramVec"]]$miniTickFactorVec[as.numeric(sel)]))
    updateNumericInput(session, "nsmall", value = as.numeric(values[["paramVec"]]$nsmallVec[as.numeric(sel)]))

    updateCheckboxInput(session, "chrSizeMbp", value = as.logical((values[["paramVec"]]$chrSizeMbpVec[as.numeric(sel)])))
    updateCheckboxInput(session, "chrNameUp", value = as.logical(values[["paramVec"]]$chrNameUpVec[as.numeric(sel)]))

    updateNumericInput(session, "distTextChr", value = as.numeric((values[["paramVec"]]$distTextChrVec[as.numeric(sel)])))
    updateNumericInput(session, "OTUTextSize", value = as.numeric((values[["paramVec"]]$OTUTextSizeVec[as.numeric(sel)])))

    updateCheckboxInput(session, "chromatids", value = values[["paramVec"]]$chromatidsVec[as.numeric(sel)])

    updateCheckboxInput(session, "holocenNotAsChromatids", value = values[["paramVec"]]$holocenNotAsChromatidsVec[as.numeric(sel)])

    updateNumericInput(session, "xModifier", value = as.numeric(values[["paramVec"]]$xModifierVec[as.numeric(sel)]))

    updateCheckboxInput(session, "circularPlot", value = values[["paramVec"]]$circularPlotVec[as.numeric(sel)])

    updateCheckboxInput(session, "ruler", value = values[["paramVec"]]$rulerVec[as.numeric(sel)])
    updateNumericInput(session, "rulerPos", value = as.numeric((values[["paramVec"]]$rulerPosVec[as.numeric(sel)])))
    updateNumericInput(session, "threshold", value = as.numeric((values[["paramVec"]]$thresholdVec[as.numeric(sel)])))

    updateNumericInput(session, "ceilingFactor", value = as.numeric((values[["paramVec"]]$ceilingFactorVec[as.numeric(sel)])))

    updateNumericInput(session, "rulerInterval", value = as.numeric((values[["paramVec"]]$rulerIntervalVec[as.numeric(sel)])))
    updateNumericInput(session, "rulerIntervalMb", value = as.numeric((values[["paramVec"]]$rulerIntervalMbVec[as.numeric(sel)])))
    updateNumericInput(session, "rulerIntervalcM", value = as.numeric((values[["paramVec"]]$rulerIntervalcMVec[as.numeric(sel)])))

    updateNumericInput(session, "ruler.tck", value = as.numeric((values[["paramVec"]]$ruler.tckVec[as.numeric(sel)])))
    updateCheckboxInput(session, "collapseCen", value = as.logical(values[["paramVec"]]$collapseCenVec[as.numeric(sel)]))

    updateNumericInput(session, "rulerNumberSize", value = as.numeric((values[["paramVec"]]$rulerNumberSizeVec[as.numeric(sel)])))

    updateNumericInput(session, "rulerNumberPos", value = as.numeric((values[["paramVec"]]$rulerNumberPosVec[as.numeric(sel)])))

    updateNumericInput(session, "rulerTitleSize", value = values[["paramVec"]]$rulerTitleSizeVec[as.numeric(sel)])

    updateNumericInput(session, "xPosRulerTitle", value = as.numeric((values[["paramVec"]]$xPosRulerTitleVec[as.numeric(sel)])))

    updateRadioButtons(session, "legend", selected = values[["paramVec"]]$legendVec[as.numeric(sel)])
    updateRadioButtons(session, "cenFormat", selected = values[["paramVec"]]$cenFormatVec[as.numeric(sel)])
    updateNumericInput(session, "cenFactor", value = values[["paramVec"]]$cenFactorVec[as.numeric(sel)])

    updateNumericInput(session, "centromereSize", value = values[["paramVec"]]$centromereSizeVec[as.numeric(sel)])
    updateCheckboxInput(session, "autoCenSize", value = as.logical(values[["paramVec"]]$autoCenSizeVec[as.numeric(sel)]))

    updateNumericInput(session, "legendWidth", value = as.numeric((values[["paramVec"]]$legendWidthVec[as.numeric(sel)])))
    updateNumericInput(session, "legendHeight", value = as.numeric((values[["paramVec"]]$legendHeightVec[as.numeric(sel)])))
    updateTextInput(session, "pattern", value = values[["paramVec"]]$patternVec[as.numeric(sel)])
    updateTextInput(session, "markNewLine", value = values[["paramVec"]]$markNewLineVec[as.numeric(sel)])
    updateTextInput(session, "forbiddenMark", value = values[["paramVec"]]$forbiddenMarkVec[as.numeric(sel)])

    updateTextInput(session, "classGroupName", value = values[["paramVec"]]$classGroupNameVec[as.numeric(sel)])

    updateCheckboxInput(session, "remSimiMarkLeg", value = as.logical((values[["paramVec"]]$remSimiMarkLegVec[as.numeric(sel)])))

    updateNumericInput(session, "markLabelSize", value = as.numeric((values[["paramVec"]]$markLabelSizeVec[as.numeric(sel)])))
    updateNumericInput(session, "markLabelSpacer", value = as.numeric((values[["paramVec"]]$markLabelSpacerVec[as.numeric(sel)])))
    updateNumericInput(session, "legendYcoord", value = as.numeric((values[["paramVec"]]$legendYcoordVec[as.numeric(sel)])))

    updateCheckboxInput(session, "fixCenBorder", value = as.logical((values[["paramVec"]]$fixCenBorderVec[as.numeric(sel)])))

    updateCheckboxInput(session, "bMarkNameAside", value = as.logical(values[["paramVec"]]$bMarkNameAsideVec[as.numeric(sel)]))

    updateCheckboxInput(session, "chrBorderColor", value = values[["paramVec"]]$chrBorderColorVec[as.numeric(sel)])
    updateNumericInput(session, "lwd.chr", value = as.numeric((values[["paramVec"]]$lwd.chrVec[as.numeric(sel)])))
    updateNumericInput(session, "lwd.cM", value = as.numeric((values[["paramVec"]]$lwd.cMVec[as.numeric(sel)])))
    updateNumericInput(session, "lwd.mimicCen", value = as.numeric((values[["paramVec"]]$lwd.mimicCenVec[as.numeric(sel)])))

    updateNumericInput(session, "xlimLeftMod", value = as.numeric((values[["paramVec"]]$xlimLeftModVec[as.numeric(sel)])))
    updateNumericInput(session, "xlimRightMod", value = values[["paramVec"]]$xlimRightModVec[as.numeric(sel)])
    updateNumericInput(session, "ylimBotMod", value = values[["paramVec"]]$ylimBotModVec[as.numeric(sel)])
    updateNumericInput(session, "ylimTopMod", value = as.numeric((values[["paramVec"]]$ylimTopModVec[as.numeric(sel)])))

    updateNumericInput(session, "hwModifier", value = as.numeric(values[["paramVec"]]$hwModifierVec[as.numeric(sel)]))

    updateNumericInput(session, "widFactor", value = as.numeric((values[["paramVec"]]$widFactorVec[as.numeric(sel)])))
    updateNumericInput(session, "heiFactor", value = as.numeric((values[["paramVec"]]$heiFactorVec[as.numeric(sel)])))

    updateRadioButtons(session, "pngorsvg", selected = values[["paramVec"]]$pngorsvgVec[as.numeric(sel)])
    updateRadioButtons(session, "pngorsvgDown", selected = values[["paramVec"]]$pngorsvgDownVec[as.numeric(sel)])

    updateTextInput(session, "mycolors", value = values[["paramVec"]]$mycolorsVec[as.numeric(sel)])
    updateTextInput(session, "markPer", value = values[["paramVec"]]$markPerVec[as.numeric(sel)])
    updateTextInput(session, "bToRemove", value = values[["paramVec"]]$bToRemoveVec[as.numeric(sel)])

    updateCheckboxInput(session, "perAsFraction", value = values[["paramVec"]]$perAsFractionVec[as.numeric(sel)])

    updateTextInput(session, "chrNamesToSwap", value = values[["paramVec"]]$chrNamesToSwapVec[as.numeric(sel)])

    updateCheckboxInput(session, "addOTUName", value = values[["paramVec"]]$addOTUNameVec[as.numeric(sel)])
    updateRadioButtons(session, "OTUfont", selected = as.character(values[["paramVec"]]$OTUfontVec[as.numeric(sel)]))

    updateTextInput(session, "moveKarHor", value = values[["paramVec"]]$moveKarHorVec[as.numeric(sel)])
    updateNumericInput(session, "mkhValue", value = values[["paramVec"]]$mkhValueVec[as.numeric(sel)])
    updateCheckboxInput(session, "anchor", value = values[["paramVec"]]$anchorVec[as.numeric(sel)])
    updateNumericInput(session, "moveAnchorV", value = values[["paramVec"]]$moveAnchorVVec[as.numeric(sel)])
    updateNumericInput(session, "moveAnchorH", value = values[["paramVec"]]$moveAnchorHVec[as.numeric(sel)])
    updateNumericInput(session, "anchorVsizeF", value = values[["paramVec"]]$anchorVsizeFVec[as.numeric(sel)])
    updateNumericInput(session, "anchorHsizeF", value = values[["paramVec"]]$anchorHsizeFVec[as.numeric(sel)])

    updateNumericInput(session, "notesTextSize", value = values[["paramVec"]]$notesTextSizeVec[as.numeric(sel)])
    updateNumericInput(session, "notesPosX", value = values[["paramVec"]]$notesPosXVec[as.numeric(sel)])

    updateNumericInput(session, "leftNoteFontUp", value = values[["paramVec"]]$leftNoteFontUpVec[as.numeric(sel)])
    updateNumericInput(session, "leftNotesPosX", value = values[["paramVec"]]$leftNotesPosXVec[as.numeric(sel)])
    updateNumericInput(session, "leftNotesPosY", value = values[["paramVec"]]$leftNotesPosYVec[as.numeric(sel)])
    updateNumericInput(session, "leftNotesUpPosY", value = values[["paramVec"]]$leftNotesUpPosYVec[as.numeric(sel)])

    updateNumericInput(session, "leftNotesUpPosX", value = values[["paramVec"]]$leftNotesUpPosXVec[as.numeric(sel)])
    updateNumericInput(session, "notesPosY", value = values[["paramVec"]]$notesPosYVec[as.numeric(sel)])

    updateNumericInput(session, "leftNoteFont", value = values[["paramVec"]]$leftNoteFontVec[as.numeric(sel)])
    updateNumericInput(session, "noteFont", value = values[["paramVec"]]$noteFontVec[as.numeric(sel)])

    updateNumericInput(session, "leftNotesUpTextSize", value = values[["paramVec"]]$leftNotesUpTextSizeVec[as.numeric(sel)])
    updateNumericInput(session, "leftNotesTextSize", value = values[["paramVec"]]$leftNotesTextSizeVec[as.numeric(sel)])

    updateCheckboxInput(session, "parseStr2lang", value = values[["paramVec"]]$parseStr2langVec[as.numeric(sel)])
    updateNumericInput(session, "moveAllKarValueHor", value = values[["paramVec"]]$moveAllKarValueHorVec[as.numeric(sel)])
    updateNumericInput(session, "moveAllKarValueY", value = values[["paramVec"]]$moveAllKarValueYVec[as.numeric(sel)])

    updateCheckboxInput(session, "verticalPlot", value = values[["paramVec"]]$verticalPlotVec[as.numeric(sel)])
    updateNumericInput(session, "karSpaceHor", value = values[["paramVec"]]$karSpaceHorVec[as.numeric(sel)])
    updateTextInput(session, "karAnchorLeft", value = values[["paramVec"]]$karAnchorLeftVec[as.numeric(sel)])

    updateCheckboxInput(session, "OTUasNote", value = values[["paramVec"]]$OTUasNoteVec[as.numeric(sel)])
    updateCheckboxInput(session, "labelOutwards", value = values[["paramVec"]]$labelOutwardsVec[as.numeric(sel)])

    # updateCheckboxInput(session, "callPlot", value = values[["paramVec"]]$callPlotVec[as.numeric(sel)] )

    updateTextInput(session, "colorBorderMark", value = values[["paramVec"]]$colorBorderMarkVec[as.numeric(sel)])
    updateNumericInput(session, "lwd.marks", value = values[["paramVec"]]$lwd.marksVec[as.numeric(sel)])
    updateCheckboxInput(session, "gishCenBorder", value = values[["paramVec"]]$gishCenBorderVec[as.numeric(sel)])
    updateNumericInput(session, "hideCenLines", value = values[["paramVec"]]$hideCenLinesVec[as.numeric(sel)])

    #
    #   nuc marks
    #
    # this one is a list #
    updateCheckboxGroupInput(session, "markType", selected = values[["paramVec"]]$markTypeVec[[as.numeric(sel)]])

    updateCheckboxGroupInput(session, "fetchSelect", selected = values[["paramVec"]]$fetchSelectVec[as.numeric(sel)])
    updateNumericInput(session, "amountofSpaces", value = values[["paramVec"]]$amountofSpacesVec[as.numeric(sel)])
    updateNumericInput(session, "colNumber", value = values[["paramVec"]]$colNumberVec[as.numeric(sel)])
    updateNumericInput(session, "protrudingInt", value = values[["paramVec"]]$protrudingIntVec[as.numeric(sel)])
    updateTextInput(session, "mycolors2", value = values[["paramVec"]]$mycolors2Vec[as.numeric(sel)])
    updateTextInput(session, "term", value = values[["paramVec"]]$termVec[as.numeric(sel)])

    updateCheckboxInput(session, "useGeneNames", value = values[["paramVec"]]$useGeneNamesVec[as.numeric(sel)])

    updateCheckboxInput(session, "useRCNames", value = values[["paramVec"]]$useRCNamesVec[as.numeric(sel)])

    updateCheckboxInput(session, "makeUnique", value = values[["paramVec"]]$makeUniqueVec[as.numeric(sel)])

    updateCheckboxInput(session, "colorFeature", value = values[["paramVec"]]$colorFeatureVec[as.numeric(sel)])

    updateRadioButtons(session, "nucMarkStyle", selected = values[["paramVec"]]$nucMarkStyleVec[as.numeric(sel)])

    updateRadioButtons(session, "pseudo", selected = values[["paramVec"]]$pseudoVec[as.numeric(sel)])

    updateCheckboxInput(session, "mirror", value = values[["paramVec"]]$mirrorVec[as.numeric(sel)])

    updateCheckboxInput(session, "addSTARTPos", value = values[["paramVec"]]$addSTARTPosVec[as.numeric(sel)])
  }
)
