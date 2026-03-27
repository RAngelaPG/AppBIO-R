# output download dataframes

#
#   download posCalc
#

output$posCalcbuttontable <- renderUI({
  req(nrow(values[["posCalc"]]) > 0) # was df
  downloadButton("downloadCsvposCalc", "Download table as .csv")
})

output$clipposCalc <- renderUI({
  req(nrow(values[["posCalc"]]) > 0) # was df
  actionButton("clipposCalcbtn", "Copy table", icon("clipboard"))
})

observeEvent(input$clipposCalcbtn, {
  clipr::write_clip(values[["posCalc"]])
})

#
#   download perMark
#

output$perMarkbuttontable <- renderUI({
  req(nrow(values[["perMark"]]) > 0) # was df
  downloadButton("downloadCsvperMark", "Download table as .csv")
})

output$clipperMark <- renderUI({
  req(nrow(values[["perMark"]]) > 0) # was df
  actionButton("clipperMarkbtn", "Copy table", icon("clipboard"))
})

observeEvent(input$clipperMarkbtn, {
  clipr::write_clip(values[["perMark"]])
})


#
#   download ARCI
#

output$ARCIbuttontable <- renderUI({
  req(nrow(values[["ARCI"]]) > 0) # was df
  downloadButton("downloadCsvARCI", "Download table as .csv")
})

output$clipARCI <- renderUI({
  req(nrow(values[["ARCI"]]) > 0) # was df
  actionButton("clipARCIbtn", "Copy table", icon("clipboard"))
})

observeEvent(input$clipARCIbtn, {
  clipr::write_clip(values[["ARCI"]])
})

#
#   download AA2
#

output$AA2buttontable <- renderUI({
  req(nrow(values[["AA2"]]) > 0) # was df
  downloadButton("downloadCsvAA2", "Download table as .csv")
})

output$clipAA2 <- renderUI({
  req(nrow(values[["AA2"]]) > 0) # was df
  actionButton("clipAA2btn", "Copy table", icon("clipboard"))
})

observeEvent(input$clipAA2btn, {
  clipr::write_clip(values[["AA2"]])
})

#
#   download chr.
#

output$buttontable <- renderUI({
  validate(
    need(try(values[["df1"]]), "") # was df
  )
  downloadButton("downloadCsv", "Download table as .csv")
})

output$buttontable2 <- renderUI({
  validate(
    need(try(values[["df1"]]), "")
  )
  if (input$saveType == "csv") {
    downloadButton("downloadCsv2", "Download table as .csv")
  } else {
    downloadButton("downloadRds2", "Download table as .rds")
  }
})

output$buttonList <- renderUI({
  validate(
    need(try(values[["dfList"]]), "")
  )
  if (input$saveType == "csv") {
    downloadButton("downloadCsv2", "Download table as .csv")
  } else {
    downloadButton("downloadRds2", "Download table as .rds")
  }
})

#
#   download notes
#

output$buttontablenotes <- renderUI({
  validate(
    need(try(values[["notes"]]), "") # was df
  )
  downloadButton("downloadCsvNotes", "Download table as .csv")
})

output$buttontablenotes2 <- renderUI({
  validate(
    need(try(values[["notes"]]), "") # wdf
  )
  if (input$saveType == "csv") {
    downloadButton("downloadCsvNotes2", "Download table as .csv")
  } else {
    downloadButton("downloadRdsNotes2", "Download table as .rds")
  }
})

#
#   download leftNotes
#

output$buttontableleftNotes <- renderUI({
  validate(
    need(try(values[["leftNotes"]]), "") # was df
  )
  downloadButton("downloadCsvleftNotes", "Download table as .csv")
})

output$buttontableleftNotes2 <- renderUI({
  validate(
    need(try(values[["leftNotes"]]), "") # wdf
  )
  if (input$saveType == "csv") {
    downloadButton("downloadCsvleftNotes2", "Download table as .csv")
  } else {
    downloadButton("downloadRdsleftNotes2", "Download table as .rds")
  }
})

#
#   download leftNotesUp
#

output$buttontableleftNotesUp <- renderUI({
  validate(
    need(try(values[["leftNotesUp"]]), "") # was df
  )
  downloadButton("downloadCsvleftNotesUp", "Download table as .csv")
})

output$buttontableleftNotesUp2 <- renderUI({
  validate(
    need(try(values[["leftNotesUp"]]), "") # wdf
  )
  if (input$saveType == "csv") {
    downloadButton("downloadCsvleftNotesUp2", "Download table as .csv")
  } else {
    downloadButton("downloadRdsleftNotesUp2", "Download table as .rds")
  }
})

#
#   download mark
#

output$buttontableMark <- renderUI({
  validate(
    need(try(values[["df1Mark"]]), "")
  )
  if (nrow(values[["df1Mark"]]) > 0) {
    downloadButton("downloadCsvMark", "Download table as .csv")
  }
})

output$buttontableMark2 <- renderUI({
  validate(
    need(try(values[["df1Mark"]]), "")
  )
  if (nrow(values[["df1Mark"]]) > 0) {
    if (input$saveType == "csv") {
      downloadButton("downloadCsvMark2", "Download table as .csv")
    } else {
      downloadButton("downloadRdsMark2", "Download table as .rds")
    }
  }
})

#
#   download mark style
#

output$buttontableMStyle <- renderUI({
  validate(
    need(try(values[["df1MStyle"]]), "")
  )
  if (nrow(values[["df1MStyle"]]) > 0) {
    downloadButton("downloadCsvMStyle", "Download table as .csv")
  }
})

output$buttontableMStyle2 <- renderUI({
  validate(
    need(try(values[["df1MStyle"]]), "")
  )
  if (nrow(values[["df1MStyle"]]) > 0) {
    if (input$saveType == "csv") {
      downloadButton("downloadCsvMStyle2", "Download table as .csv")
    } else {
      downloadButton("downloadRdsMStyle2", "Download table as .rds")
    }
  }
})

#
#   downloadHandler perMark
#

output$downloadCsvperMark <- downloadHandler(
  filename = function() {
    # paste0(input$chrfilename)
    paste0(input$perMarkfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["perMark"]], file,
      na = "", # wdf
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

#
#   downloadHandler posCalc
#

output$downloadCsvposCalc <- downloadHandler(
  filename = function() {
    # paste0(input$chrfilename)
    paste0(input$posCalcfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["posCalc"]], file,
      na = "", # wdf
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

#
#   downloadHandler ARCI
#

output$downloadCsvARCI <- downloadHandler(
  filename = function() {
    # paste0(input$chrfilename)
    paste0(input$ARCIfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["ARCI"]], file,
      na = "", # wdf
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

#
#   downloadHandler AA2
#

output$downloadCsvAA2 <- downloadHandler(
  filename = function() {
    # paste0(input$chrfilename)
    paste0(input$AA2filename, ".csv")
  },
  content = function(file) {
    write.csv(values[["AA2"]], file,
      na = "", # wdf
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)


#
#   downloadHandler chr
#

output$downloadCsv <- downloadHandler(
  filename = function() {
    # paste0(input$chrfilename)
    paste0(input$chrfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["df1"]], file,
      na = "", # wdf
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadCsv2 <- downloadHandler(
  filename = function() {
    paste0(input$chrfilename2, ".csv")
  },
  content = function(file) {
    write.csv(values[["df1"]], file,
      na = "", # wdf
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadRds2 <- downloadHandler(
  filename = function() {
    paste0(input$chrfilename2, ".rds")
  },
  content = function(file) {
    saveRDS(values[["df1"]], file, )
  },
  contentType = "rds"
)

#
#   downloadHandler Mark
#

output$downloadCsvMark <- downloadHandler(
  filename = function() {
    # paste0(input$markfilename)
    paste0(input$markfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["df1Mark"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadCsvMark2 <- downloadHandler(
  filename = function() {
    paste0(input$markfilename2, ".csv")
  },
  content = function(file) {
    write.csv(values[["df1Mark"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadRdsMark2 <- downloadHandler(
  filename = function() {
    paste0(input$markfilename2, ".rds")
  },
  content = function(file) {
    saveRDS(values[["df1Mark"]], file, )
  },
  contentType = "rds"
)


#
#   downloadHandler MStyle style
#

output$downloadCsvMStyle <- downloadHandler(
  filename = function() {
    # paste0(input$MStylefilename)
    paste0(input$MStylefilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["df1MStyle"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadCsvMStyle2 <- downloadHandler(
  filename = function() {
    paste0(input$MStylefilename2, ".csv")
  },
  content = function(file) {
    write.csv(values[["df1MStyle"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadRdsMStyle2 <- downloadHandler(
  filename = function() {
    paste0(input$MStylefilename2, ".rds")
  },
  content = function(file) {
    saveRDS(values[["df1MStyle"]], file, )
  },
  contentType = "rds"
)

#
#   downloadHandler notes
#

output$downloadCsvNotes <- downloadHandler(
  filename = function() {
    paste0(input$notesfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["notes"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadCsvNotes2 <- downloadHandler(
  filename = function() {
    paste0(input$notesfilename2, ".csv")
  },
  content = function(file) {
    write.csv(values[["notes"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadRdsNotes2 <- downloadHandler(
  filename = function() {
    paste0(input$notesfilename2, ".rds")
  },
  content = function(file) {
    saveRDS(values[["notes"]], file, )
  },
  contentType = "rds"
)

#
#   downloadHandler leftNotes
#

output$downloadCsvleftNotes <- downloadHandler(
  filename = function() {
    paste0(input$leftNotesfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["leftNotes"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadCsvleftNotes2 <- downloadHandler(
  filename = function() {
    paste0(input$leftNotesfilename2, ".csv")
  },
  content = function(file) {
    write.csv(values[["leftNotes"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadRdsleftNotes2 <- downloadHandler(
  filename = function() {
    paste0(input$leftNotesfilename2, ".rds")
  },
  content = function(file) {
    saveRDS(values[["leftNotes"]], file, )
  },
  contentType = "rds"
)

#
#   downloadHandler leftNotesUp
#

output$downloadCsvleftNotesUp <- downloadHandler(
  filename = function() {
    paste0(input$leftNotesUpfilename, ".csv")
  },
  content = function(file) {
    write.csv(values[["leftNotesUp"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadCsvleftNotesUp2 <- downloadHandler(
  filename = function() {
    paste0(input$leftNotesUpfilename2, ".csv")
  },
  content = function(file) {
    write.csv(values[["leftNotesUp"]], file,
      na = "",
      row.names = FALSE, quote = TRUE
    )
  },
  contentType = "text/csv"
)

output$downloadRdsleftNotesUp2 <- downloadHandler(
  filename = function() {
    paste0(input$leftNotesUpfilename2, ".rds")
  },
  content = function(file) {
    saveRDS(values[["leftNotesUp"]], file, )
  },
  contentType = "rds"
)

#
#   render posCalc table
#

output$outposCalc <- renderRHandsontable({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "marksTab")
  req(nrow(values[["posCalc"]]) > 0)
  rhandsontable(values[["posCalc"]], readOnly = TRUE) %>%
    hot_cols(columnSorting = TRUE)
})

output$outposCalcUI <- renderUI({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "marksTab")
  # wellPanel(
  rHandsontableOutput("outposCalc")
  # )
})

#
#   render perMark table
#

output$outperMark <- renderRHandsontable({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "marksTab")
  req(nrow(values[["perMark"]]) > 0)
  rhandsontable(values[["perMark"]], readOnly = TRUE) %>%
    hot_cols(columnSorting = TRUE)
})

output$outperMarkUI <- renderUI({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "marksTab")
  # wellPanel(
  rHandsontableOutput("outperMark")
  # )
})

#
#   render ARCI table
#

output$outARCI <- renderRHandsontable({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "indicesTab")
  req(nrow(values[["ARCI"]]) > 0)
  rhandsontable(values[["ARCI"]], readOnly = TRUE) %>%
    hot_cols(columnSorting = TRUE)
})

output$outARCIUI <- renderUI({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "indicesTab")
  # wellPanel(
  rHandsontableOutput("outARCI")
  # )
})

#
#   render AA2 table
#

output$outAA2 <- renderRHandsontable({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "indicesTab")
  req(nrow(values[["AA2"]]) > 0)
  rhandsontable(values[["AA2"]], readOnly = TRUE) %>%
    hot_cols(columnSorting = TRUE)
})

output$outAA2UI <- renderUI({
  req(CurrentM$menu == "indicesMenu", Current$Tab5 == "indicesTab")
  # wellPanel(
  rHandsontableOutput("outAA2")
  # )
})

#
#   render notes table
#

output$outNotes <- renderRHandsontable({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "notesTab")
  req(nrow(values[["notes"]]) > 0)
  rhandsontable(values[["notes"]]) %>%
    hot_cols(columnSorting = TRUE)
})

#
#   render leftNotes table
#

output$outleftNotes <- renderRHandsontable({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "notesTab")
  req(nrow(values[["leftNotes"]]) > 0)
  rhandsontable(values[["leftNotes"]]) %>%
    hot_cols(columnSorting = TRUE)
})

#
#   render leftNotesUp table
#

output$outleftNotesUp <- renderRHandsontable({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "notesTab")
  req(nrow(values[["leftNotesUp"]]) > 0)

  rhandsontable(values[["leftNotesUp"]]) %>%
    hot_cols(columnSorting = TRUE)
})


#
#   render chr table
#

output$out <- renderRHandsontable({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "dfChrTab")
  req(nrow(values[["df1"]]) > 0)
  rhandsontable(values[["df1"]]) %>%
    hot_cols(columnSorting = TRUE)
})

#
#   render mark table
#

output$outMark <- renderRHandsontable({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "dfMarkTab")
  req(nrow(values[["df1Mark"]]) > 0)

  rhandsontable(values[["df1Mark"]]) %>%
    hot_cols(columnSorting = TRUE)
})

#
#   render mark style table
#

output$outMStyle <- renderRHandsontable({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "dfMSTab")
  req(nrow(values[["df1MStyle"]]) > 0)

  validate(need(
    try(values[["df1MStyle"]]),
    "Hint: Colors can be added in Parameters Menu under 'Color' (optional)"
  ))
  if (invalid(values[["df1MStyle"]])) {
    # rhandsontable(values[["df1MStyle"]])
  } else if (nrow(values[["df1MStyle"]]) > 0 & "style" %in% colnames(values[["df1MStyle"]])) {
    Options <- c(
      "dots", "square", "squareLeft", "cM", "cMLeft", "cenStyle", "upArrow", "downArrow",
      "exProtein", "inProtein"
    )

    tmpExampleTable <- rhandsontable(values[["df1MStyle"]],
      rowHeaders = NULL,
      stretchH = "all",
      selectCallback = TRUE,
      width = 400, height = 400
    ) %>%
      hot_col("style",
        allowInvalid = FALSE, type = "dropdown",
        readOnly = TRUE
      )
    if (!is.null(input$outMStyle_select$select$r)) {
      tmpExampleTable <- hot_col(tmpExampleTable,
        col = "style",
        allowInvalid = FALSE,
        type = "dropdown",
        source = Options
      ) %>%
        hot_cell(
          row = input$outMStyle_select$select$r,
          col = "style",
          readOnly = FALSE
        )
    }
    tmpExampleTable %>% hot_cols(columnSorting = TRUE)
  } else if (nrow(values[["df1MStyle"]]) > 0) {
    rhandsontable(values[["df1MStyle"]])
  }
})

output$outUI <- renderUI({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "dfChrTab")
  rHandsontableOutput("out")
})

output$outnotesUI <- renderUI({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "notesTab")
  rHandsontableOutput("outNotes")
})

output$outleftNotesUI <- renderUI({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "notesTab")
  rHandsontableOutput("outleftNotes")
})

output$outleftNotesUpUI <- renderUI({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "notesTab")
  rHandsontableOutput("outleftNotesUp")
})

output$outMarkUI <- renderUI({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "dfMarkTab")
  rHandsontableOutput("outMark")
})

output$outMStyleUI <- renderUI({
  req(CurrentM$menu == "DFsMenu", Current$Tab == "dfMSTab")
  rHandsontableOutput("outMStyle")
})
