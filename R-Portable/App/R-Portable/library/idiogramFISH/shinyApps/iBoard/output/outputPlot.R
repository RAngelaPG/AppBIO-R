output$idiogramPlot <- renderImage(
  {
    validate(need(try(values[["stop"]] == FALSE), "not ready"))
    validate(need(try(values[["df1"]]), "Start with Examples, data.frames or Nucleotide pages"))
    validate(need(try(values[["go"]] == TRUE), "not ready"))
    validate(need(try(values[["df1"]]), helpString))

    validate(need(
      try(inherits(values[["df1"]], "data.frame")),
      helpString
    ))
    validate(need(
      try(nrow(values[["df1"]]) > 0),
      helpString
    ))

    Sys.sleep(0.2)

    showModal(modalDialog(
      title = "Plotting",
      tagList(
        paste("this pop-up will close after completion. Press ESC to wait in shiny app."),
        br(),
        paste0("display: ", values[["pngorsvg"]], values[["decision"]])
      ),
      easyClose = TRUE,
      footer = tagList(
        modalButton("Wait in shiny app")
      )
    ))

    info <- getCurrentOutputInfo()
    width <- as.numeric(isolate(info$width())) * (as.numeric(input$widFactor) / 100)
    height <- as.numeric(isolate(info$width())) *
      (as.numeric(input$heiFactor) * (as.numeric(input$widFactor) / 100))

    mysvgwidth <- width / 80
    mysvgheight <- height / 80

    values[["width"]] <- width
    values[["height"]] <- height
    values[["mysvgwidth"]] <- mysvgwidth
    values[["mysvgheight"]] <- mysvgheight


    if (values[["pngorsvg"]] == "svg") {
      outfileSvg <- tempfile(fileext = ".svg")

      values[["imageType"]] <- "image/svg+xml"

      svg(outfileSvg,
        width = mysvgwidth,
        height = mysvgheight
      )
    } else {
      outfilePng <- tempfile(fileext = ".png")
      values[["imageType"]] <- "image/png"
      png(outfilePng, width = width, height = height)
    }

    tryCatch(
      capture.output(
        plotIdiograms(
          dfChrSize = values[["df1"]],
          dfMarkColor = values[["df1MStyle"]],
          dfMarkPos = values[["df1Mark"]],
          notes = values[["notes"]],
          leftNotesUp = values[["leftNotesUp"]],
          leftNotes = values[["leftNotes"]],
          colorBorderMark = input$colorBorderMark,
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
          karHeight = input$karHeight,
          karHeiSpace = input$karHeiSpace,
          amoSepar = input$amoSepar,
          karSepar = input$karSepar,
          legend = as.character(input$legend),
          mycolors = if (length(values[["mycolors"]]) == 0) {
            ""
          } else {
            values[["mycolors"]]
          },
          markPer = if (length(values[["markPer"]]) == 0) {
            ""
          } else {
            values[["markPer"]]
          },
          bToRemove = if (length(values[["bToRemove"]]) == 0) {
            ""
          } else {
            values[["bToRemove"]]
          },
          bannedMarkName = if (length(values[["bannedMarkName"]]) == 0) {
            ""
          } else {
            values[["bannedMarkName"]]
          },
          specialOTUNames = if (length(values[["specialOTUNames"]]) == 0) {
            ""
          } else {
            values[["specialOTUNames"]]
          },
          missOTUspacings = if (length(values[["missOTUspacings"]]) == 0) {
            ""
          } else {
            values[["missOTUspacings"]]
          },
          addMissingOTUAfter = if (length(values[["addMissingOTUAfter"]]) == 0) {
            ""
          } else {
            values[["addMissingOTUAfter"]]
          },
          specialyTitle = input$specialyTitle,
          specialChrWidth = input$specialChrWidth,
          specialChrSpacing = input$specialChrSpacing,
          origin = input$origin,
          OTUfamily = input$OTUfamily,
          classChrName = input$classChrName,
          classChrNameUp = input$classChrNameUp,
          chrIdPatternRem = input$chrIdPatternRem,
          xModMonoHoloRate = input$xModMonoHoloRate,
          yTitle = input$yTitle,
          nameChrIndexPos = input$nameChrIndexPos,
          perAsFraction = input$perAsFraction,
          cenFormat = input$cenFormat,
          cenFactor = input$cenFactor,
          autoCenSize = input$autoCenSize,
          centromereSize = input$centromereSize,
          chrColor = input$chrColor,
          cenColor = input$cenColor,
          chrId = input$chrId,
          chrWidth = input$chrWidth,
          chrSpacing = input$chrSpacing,
          groupSepar = input$groupSepar,
          squareness = input$squareness,
          markDistType = input$markDistType,
          markN = input$markN,
          n = input$n,
          orderChr = input$orderChr,
          useOneDot = input$useOneDot,
          pMarkFac = input$pMarkFac,
          cMBeginCenter = input$cMBeginCenter,
          protruding = input$protruding,
          arrowhead = input$arrowhead,
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
          OTUplacing = input$OTUplacing,
          morpho = input$morpho,
          chrIndex = input$chrIndex,
          karIndex = input$karIndex,
          karIndexPos = input$karIndexPos,
          chrSize = input$chrSize,
          miniTickFactor = input$miniTickFactor,
          nsmall = input$nsmall,
          showMarkPos = input$showMarkPos,
          indexIdTextSize = input$indexIdTextSize,
          chrSizeMbp = input$chrSizeMbp,
          chrNameUp = input$chrNameUp,
          useMinorTicks = input$useMinorTicks,
          chromatids = input$chromatids,
          holocenNotAsChromatids = input$holocenNotAsChromatids,
          xModifier = input$xModifier,
          circularPlot = input$circularPlot,
          addOTUName = input$addOTUName,
          OTUfont = as.numeric(input$OTUfont),
          moveKarHor = as.character(input$moveKarHor),
          mkhValue = input$mkhValue,
          anchor = input$anchor,
          moveAnchorV = input$moveAnchorV,
          moveAnchorH = input$moveAnchorH,
          anchorVsizeF = input$anchorVsizeF,
          anchorHsizeF = input$anchorHsizeF,
          parseStr2lang = input$parseStr2lang,
          moveAllKarValueHor = input$moveAllKarValueHor,
          moveAllKarValueY = input$moveAllKarValueY,
          verticalPlot = input$verticalPlot,
          karSpaceHor = input$karSpaceHor,
          karAnchorLeft = input$karAnchorLeft,
          OTUasNote = input$OTUasNote,
          labelOutwards = input$labelOutwards,
          notesTextSize = input$notesTextSize,
          notesPosX = input$notesPosX,
          ruler = input$ruler,
          rulerPos = input$rulerPos,
          rulerInterval = input$rulerInterval,
          rulerIntervalMb = input$rulerIntervalMb,
          rulerIntervalcM = input$rulerIntervalcM,
          threshold = input$threshold,
          ceilingFactor = input$ceilingFactor,
          ruler.tck = input$ruler.tck,
          collapseCen = input$collapseCen,
          rulerNumberSize = input$rulerNumberSize,
          rulerNumberPos = input$rulerNumberPos,
          rulerTitleSize = input$rulerTitleSize,
          xPosRulerTitle = input$xPosRulerTitle,
          legendWidth = input$legendWidth,
          legendHeight = input$legendHeight,
          pattern = input$pattern,
          forbiddenMark = input$forbiddenMark,
          markNewLine = input$markNewLine,
          classGroupName = input$classGroupName,
          remSimiMarkLeg = input$remSimiMarkLeg,
          markLabelSize = input$markLabelSize,
          markLabelSpacer = input$markLabelSpacer,
          legendYcoord = input$legendYcoord,
          fixCenBorder = input$fixCenBorder,
          bMarkNameAside = input$bMarkNameAside,
          chrBorderColor = input$chrBorderColor,
          lwd.chr = input$lwd.chr,
          lwd.cM = input$lwd.cM,
          lwd.mimicCen = input$lwd.mimicCen,
          distTextChr = input$distTextChr,
          OTUTextSize = input$OTUTextSize,
          xlimLeftMod = input$xlimLeftMod,
          xlimRightMod = input$xlimRightMod,
          ylimBotMod = input$ylimBotMod,
          ylimTopMod = input$ylimTopMod
        ),
        file = (outfile <- file(
          filenameR(), "w"
        )), type = "message"
      ),
      error = function(e) {
        helpString
      }
    )

    close(outfile)

    dev.off()

    removeModal()

    hwModifier1 <- ifelse(input$hwModifier == 0, 0.1, input$hwModifier)

    if (values[["pngorsvg"]] == "svg") {
      list(
        src = normalizePath(outfileSvg),
        contentType = "image/svg+xml",
        width = width * hwModifier1,
        height = height * hwModifier1,
        alt = "My plot"
      )
    } else {
      list(
        src = normalizePath(outfilePng),
        contentType = "image/png",
        width = width * hwModifier1,
        height = height * hwModifier1,
        alt = "My plot"
      )
    }
  },
  deleteFile = TRUE
)

output$buttonScript <- renderUI({
  validate(
    need(try(scriptR()), "")
  )
  downloadButton("downloadR", "Download R-script")
})

output$buttonPresets <- renderUI({
  validate(
    need(try(!is.null(presetsR())), "valid plot needed")
  )
  downloadButton("downloadPresets", "Download Presets", class = "btn-success")
})

output$buttonPresets2 <- renderUI({
  validate(
    need(try(!is.null(presetsR())), "valid plot needed")
  )
  downloadButton("downloadPresets2", "Download Presets", class = "btn-success")
})

output$downloadPresets <- downloadHandler(
  filename = function() {
    "preset.rds"
  },
  content = function(file) {
    saveRDS(presetsR(), file, )

    pre <- presetsR()

    values[["upPresetName"]] <- "added"

    values[["upPresetOrigin"]] <- "file was loaded"

    current_len <- length(values[["paramVec"]][[1]])

    for (i in 1:length(values[["paramVec"]])) {
      if (length(values[["paramVec"]][[i]]) > current_len) {
        current_len <- length(values[["paramVec"]][[i]])
      }
    }

    for (i in 1:length(pre)) {
      name_pa <- sub("(.*)", "\\1Vec", attr(pre[i], "name"))
      if (length(pre[[i]]) > 1) {
        if (inherits(pre[[i]], "data.frame")) {
          values[["paramVec"]][[name_pa]][(current_len + 1)] <- list(pre[[i]])
        } else {
          values[["paramVec"]][[name_pa]][(current_len + 1)] <- paste(pre[[i]], collapse = ",")
        }
      } else {
        values[["paramVec"]][[name_pa]][(current_len + 1)] <- pre[[i]]
      }
    }

    for (i in 1:length(values[["paramVec"]])) {
      if (length(values[["paramVec"]][[i]]) > current_len) {
        current_len <- length(values[["paramVec"]][[i]])
      }
    }

    values[["maxEx"]] <- values[["current_len"]] <- current_len

    for (i in 1:length(values[["paramVec"]])) {
      par_name <- attr(values[["paramVec"]][i], "name")
      name_def <- sub("(.*)Vec", "\\1Default", par_name)
      # complete missing
      if (length(values[["paramVec"]][[i]]) < current_len) {
        if (length(paramValues[[name_def]]) > 1) {
          #
          # marktype is longer than 1
          #
          values[["paramVec"]][[par_name]][(current_len)] <- list(paramValues[[name_def]])
        } else {
          values[["paramVec"]][[par_name]][(current_len)] <- paramValues[[name_def]]
        }
      }
    }
    values[["presetStatus"]] <- "Ready"
  },
  contentType = "rds"
)

output$downloadPresets2 <- downloadHandler(
  filename = function() {
    "preset.rds"
  },
  content = function(file) {
    saveRDS(presetsR(), file, )

    pre <- presetsR()

    values[["upPresetName"]] <- "added"

    values[["upPresetOrigin"]] <- "file was loaded"

    current_len <- length(values[["paramVec"]][[1]])

    for (i in 1:length(values[["paramVec"]])) {
      if (length(values[["paramVec"]][[i]]) > current_len) {
        current_len <- length(values[["paramVec"]][[i]])
      }
    }

    for (i in 1:length(pre)) {
      name_pa <- sub("(.*)", "\\1Vec", attr(pre[i], "name"))
      if (length(pre[[i]]) > 1) {
        if (inherits(pre[[i]], "data.frame")) {
          values[["paramVec"]][[name_pa]][(current_len + 1)] <- list(pre[[i]])
        } else {
          values[["paramVec"]][[name_pa]][(current_len + 1)] <- paste(pre[[i]], collapse = ",")
        }
      } else {
        values[["paramVec"]][[name_pa]][(current_len + 1)] <- pre[[i]]
      }
    }

    for (i in 1:length(values[["paramVec"]])) {
      if (length(values[["paramVec"]][[i]]) > current_len) {
        current_len <- length(values[["paramVec"]][[i]])
      }
    }

    values[["maxEx"]] <- values[["current_len"]] <- current_len

    for (i in 1:length(values[["paramVec"]])) {
      par_name <- attr(values[["paramVec"]][i], "name")
      name_def <- sub("(.*)Vec", "\\1Default", par_name)
      # complete missing
      if (length(values[["paramVec"]][[i]]) < current_len) {
        if (length(paramValues[[name_def]]) > 1) {
          #
          # marktype is longer than 1
          #
          values[["paramVec"]][[par_name]][(current_len)] <- list(paramValues[[name_def]])
        } else {
          values[["paramVec"]][[par_name]][(current_len)] <- paramValues[[name_def]]
        }
      }
    }
    values[["presetStatus"]] <- "Ready"
  },
  contentType = "rds"
)

output$downloadR <- downloadHandler(
  filename = "script.R",
  content = function(file) {
    writeLines(scriptR(), file)
  }
)

output$pngorsvgDownButton <- downloadHandler(
  filename <- function() {
    paste0("plot.", input$pngorsvgDown)
  },
  content = function(file) {
    if (input$pngorsvgDown == "svg" & values[["number"]] == 10) {
      if (Sys.info()["sysname"] == "Linux") {
        message1 <- "Linux might NOT be able to display .svg with more than 200,000 polygons"
        message2 <- "Search for a workaround to increase this number for your system"
      } else if (Sys.info()["sysname"] == "Windows") {
        message1 <- ""
        message2 <- ""
      } else {
        message1 <- "Problems known in linux for opening large .svg"
        message2 <- "You may have to search for a workaround for your system"
      }

      showModal(modalDialog(
        title = "Plotting to file, please wait.",
        tagList(
          "You selected to download a large chr. in .svg format",
          br(),
          message1,
          br(),
          message2,
          br(),
          "this pop-up will close after completion. Press ESC to wait in shiny app"
        ),
        easyClose = TRUE,
        footer = list(
          modalButton("Wait in shiny app")
        )
      ))
    } else {
      showModal(modalDialog(
        title = "Plotting to file, please wait",
        "this pop-up will close after completion. Press ESC to wait in shiny app",
        easyClose = TRUE,
        footer = modalButton("Wait in shiny app")
      ))
    }

    if (input$pngorsvgDown == "svg") {
      svg(file, width = values[["mysvgwidth"]], height = values[["mysvgheight"]])
    } else {
      png(file, width = values[["width"]], height = values[["height"]])
    }

    plotIdiograms(
      dfChrSize = values[["df1"]],
      dfMarkColor = values[["df1MStyle"]],
      dfMarkPos = values[["df1Mark"]],
      notes = values[["notes"]],
      leftNotes = values[["leftNotes"]],
      leftNotesUp = values[["leftNotesUp"]],
      colorBorderMark = input$colorBorderMark,
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
      karHeight = input$karHeight,
      karHeiSpace = input$karHeiSpace,
      amoSepar = input$amoSepar,
      karSepar = input$karSepar,
      legend = as.character(input$legend),
      mycolors = if (length(values[["mycolors"]]) == 0) {
        ""
      } else {
        values[["mycolors"]]
      },
      markPer = if (length(values[["markPer"]]) == 0) {
        ""
      } else {
        values[["markPer"]]
      },
      bToRemove = if (length(values[["bToRemove"]]) == 0) {
        ""
      } else {
        values[["bToRemove"]]
      },
      bannedMarkName = if (length(values[["bannedMarkName"]]) == 0) {
        ""
      } else {
        values[["bannedMarkName"]]
      },
      specialOTUNames = if (length(values[["specialOTUNames"]]) == 0) {
        ""
      } else {
        values[["specialOTUNames"]]
      },
      missOTUspacings = if (length(values[["missOTUspacings"]]) == 0) {
        ""
      } else {
        values[["missOTUspacings"]]
      },
      addMissingOTUAfter = if (length(values[["addMissingOTUAfter"]]) == 0) {
        ""
      } else {
        values[["addMissingOTUAfter"]]
      },
      specialyTitle = input$specialyTitle,
      specialChrWidth = input$specialChrWidth,
      specialChrSpacing = input$specialChrSpacing,
      origin = input$origin,
      OTUfamily = input$OTUfamily,
      classChrName = input$classChrName,
      classChrNameUp = input$classChrNameUp,
      chrIdPatternRem = input$chrIdPatternRem,
      xModMonoHoloRate = input$xModMonoHoloRate,
      yTitle = input$yTitle,
      nameChrIndexPos = input$nameChrIndexPos,
      perAsFraction = input$perAsFraction,
      cenFormat = input$cenFormat,
      cenFactor = input$cenFactor,
      autoCenSize = input$autoCenSize,
      centromereSize = input$centromereSize,
      chrColor = input$chrColor,
      cenColor = input$cenColor,
      chrId = input$chrId,
      chrWidth = input$chrWidth,
      chrSpacing = input$chrSpacing,
      groupSepar = input$groupSepar,
      squareness = input$squareness,
      markDistType = input$markDistType,
      markN = input$markN,
      n = input$n,
      orderChr = input$orderChr,
      useOneDot = input$useOneDot,
      pMarkFac = input$pMarkFac,
      cMBeginCenter = input$cMBeginCenter,
      protruding = input$protruding,
      arrowhead = input$arrowhead,
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
      OTUplacing = input$OTUplacing,
      morpho = input$morpho,
      chrIndex = input$chrIndex,
      karIndex = input$karIndex,
      karIndexPos = input$karIndexPos,
      chrSize = input$chrSize,
      miniTickFactor = input$miniTickFactor,
      nsmall = input$nsmall,
      showMarkPos = input$showMarkPos,
      indexIdTextSize = input$indexIdTextSize,
      chrSizeMbp = input$chrSizeMbp,
      chrNameUp = input$chrNameUp,
      useMinorTicks = input$useMinorTicks,
      chromatids = input$chromatids,
      holocenNotAsChromatids = input$holocenNotAsChromatids,
      xModifier = input$xModifier,
      circularPlot = input$circularPlot,
      addOTUName = input$addOTUName,
      OTUfont = as.numeric(input$OTUfont),
      moveKarHor = as.character(input$moveKarHor),
      mkhValue = input$mkhValue,
      anchor = input$anchor,
      moveAnchorV = input$moveAnchorV,
      moveAnchorH = input$moveAnchorH,
      anchorVsizeF = input$anchorVsizeF,
      anchorHsizeF = input$anchorHsizeF,
      parseStr2lang = input$parseStr2lang,
      moveAllKarValueHor = input$moveAllKarValueHor,
      moveAllKarValueY = input$moveAllKarValueY,
      verticalPlot = input$verticalPlot,
      karSpaceHor = input$karSpaceHor,
      karAnchorLeft = input$karAnchorLeft,
      OTUasNote = input$OTUasNote,
      labelOutwards = input$labelOutwards,
      notesTextSize = input$notesTextSize,
      notesPosX = input$notesPosX,
      ruler = input$ruler,
      rulerPos = input$rulerPos,
      rulerInterval = input$rulerInterval,
      rulerIntervalMb = input$rulerIntervalMb,
      rulerIntervalcM = input$rulerIntervalcM,
      threshold = input$threshold,
      ceilingFactor = input$ceilingFactor,
      ruler.tck = input$ruler.tck,
      collapseCen = input$collapseCen,
      rulerNumberSize = input$rulerNumberSize,
      rulerNumberPos = input$rulerNumberPos,
      rulerTitleSize = input$rulerTitleSize,
      xPosRulerTitle = input$xPosRulerTitle,
      legendWidth = input$legendWidth,
      legendHeight = input$legendHeight,
      pattern = input$pattern,
      forbiddenMark = input$forbiddenMark,
      markNewLine = input$markNewLine,
      classGroupName = input$classGroupName,
      remSimiMarkLeg = input$remSimiMarkLeg,
      markLabelSize = input$markLabelSize,
      markLabelSpacer = input$markLabelSpacer,
      legendYcoord = input$legendYcoord,
      fixCenBorder = input$fixCenBorder,
      bMarkNameAside = input$bMarkNameAside,
      chrBorderColor = input$chrBorderColor,
      lwd.chr = input$lwd.chr,
      lwd.cM = input$lwd.cM,
      lwd.mimicCen = input$lwd.mimicCen,
      distTextChr = input$distTextChr,
      OTUTextSize = input$OTUTextSize,
      xlimLeftMod = input$xlimLeftMod,
      xlimRightMod = input$xlimRightMod,
      ylimBotMod = input$ylimBotMod,
      ylimTopMod = input$ylimTopMod
    )

    dev.off()
    removeModal()
  },
  contentType = values[["imageType"]]
)
