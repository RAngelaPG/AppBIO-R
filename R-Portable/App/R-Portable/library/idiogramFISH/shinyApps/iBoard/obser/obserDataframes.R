#
#   loadDFbutton is a Nuccore button
#

observeEvent(input$loadDFbutton, ignoreInit = T, {
  values[["df1MStyleName"]] <- values[["df1MarkName"]] <- values[["df1Name"]] <- "Nucleotide tab"
  values[["df1MStyleOrigin"]] <- values[["df1MarkOrigin"]] <- values[["df1Origin"]] <- "data.frame was loaded"
})

# update all origins

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0), "not ready"))
    Sys.sleep(.1)
    values[["df1MStyleOrigin"]] <- values[["df1MarkOrigin"]] <- values[["df1Origin"]] <- "data.frame was loaded"
    values[["notesOrigin"]] <- "data.frame was loaded"
    values[["leftNotesOrigin"]] <- "data.frame was loaded"
    values[["leftNotesUpOrigin"]] <- "data.frame was loaded"
  }
)

#
#   update chr. data.frame
#

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0 & numberRe() <= maxEx), "not ready"))
    # values[["go"]]<-FALSE
    if (values[["canButton"]]) {
      values[["df1"]] <- data.frame()
      Sys.sleep(0.1)

      df <- values[["paramVec"]]$dfChrSizeVec[[as.numeric(numberRe())]]
      tryCatch(load(paste0("www/rda/", df, ".rda")), error = function(e) {
        ""
      }, warning = function(w) {
        ""
      })
      values[["df1"]] <- get(df)
      values[["df1Name"]] <- ifelse(df == "emptydata.frame", "no", df)
    }
    # values[["go"]]<-TRUE
  }
)

observeEvent(input$file1, {
  inFile <- input[["file1"]]
  if (grepl(".csv", inFile$datapath)) {
    values[["df1"]] <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(".rds", inFile$datapath)) {
    values[["df1"]] <- readRDS(inFile$datapath)
  }
  values[["df1Name"]] <- inFile$name
  values[["df1Origin"]] <- "file was loaded"
})

observeEvent(input$upFilePreset, {
  inFile <- input[["upFilePreset"]]

  if (grepl(".rds", inFile$datapath)) {
    pre <- readRDS(inFile$datapath)

    values[["upPresetName"]] <- inFile$name

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
    values[["pre"]] <- pre
  } else {
    values[["upPresetOrigin"]] <- "file must be .rds"
  }
})

#
#   update notes data.frame
#

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0 & numberRe() <= maxEx), "not ready"))
    if (values[["canButton"]]) {
      values[["notes"]] <- data.frame()
      Sys.sleep(0.1)
      df <- values[["paramVec"]]$notesVec[[as.numeric(numberRe())]]
      tryCatch(load(paste0("www/rda/", df, ".rda")), error = function(e) {
        ""
      }, warning = function(w) {
        ""
      })
      values[["notes"]] <- get(df)
      values[["notesName"]] <- ifelse(df == "emptydata.frame", "no", df)
    }
  }
)

observeEvent(input$file1Notes, {
  inFile <- input[["file1Notes"]]

  if (grepl(".csv", inFile$datapath)) {
    values[["notes"]] <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(".rds", inFile$datapath)) {
    values[["notes"]] <- readRDS(inFile$datapath)
  }

  values[["notesName"]] <- inFile$name
  values[["notesOrigin"]] <- "file was loaded"
})

#
#   update leftNotes data.frame
#

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0 & numberRe() <= maxEx), "not ready"))
    if (values[["canButton"]]) {
      values[["leftNotes"]] <- data.frame()
      Sys.sleep(0.1)
      df <- values[["paramVec"]]$leftNotesVec[[as.numeric(numberRe())]]
      tryCatch(load(paste0("www/rda/", df, ".rda")), error = function(e) {
        ""
      }, warning = function(w) {
        ""
      })
      values[["leftNotes"]] <- get(df)
      values[["leftNotesName"]] <- ifelse(df == "emptydata.frame", "no", df)
    }
  }
)

observeEvent(input$file1leftNotes, {
  inFile <- input[["file1leftNotes"]]

  if (grepl(".csv", inFile$datapath)) {
    values[["leftNotes"]] <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(".rds", inFile$datapath)) {
    values[["leftNotes"]] <- readRDS(inFile$datapath)
  }

  values[["leftNotesName"]] <- inFile$name
  values[["leftNotesOrigin"]] <- "file was loaded"
})

#
#   update leftNotesUp data.frame
#

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0 & numberRe() <= maxEx), "not ready"))
    if (values[["canButton"]]) {
      values[["leftNotesUp"]] <- data.frame()
      Sys.sleep(0.1)
      df <- values[["paramVec"]]$leftNotesUpVec[[as.numeric(numberRe())]]
      tryCatch(load(paste0("www/rda/", df, ".rda")), error = function(e) {
        ""
      }, warning = function(w) {
        ""
      })
      values[["leftNotesUp"]] <- get(df)
      values[["leftNotesUpName"]] <- ifelse(df == "emptydata.frame", "no", df)
    }
  }
)

observeEvent(input$file1leftNotesUp, {
  inFile <- input[["file1leftNotesUp"]]

  if (grepl(".csv", inFile$datapath)) {
    values[["leftNotesUp"]] <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(".rds", inFile$datapath)) {
    values[["leftNotesUp"]] <- readRDS(inFile$datapath)
  }

  values[["leftNotesUpName"]] <- inFile$name
  values[["leftNotesUpOrigin"]] <- "file was loaded"
})


#
#   update mark pos data.frame
#

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0 & numberRe() <= maxEx), "not ready"))
    # par228<<-values[["paramVec"]]
    if (values[["canButton"]]) {
      values[["df1Mark"]] <- data.frame()
      Sys.sleep(0.1)
      df <- values[["paramVec"]]$dfMarkPosVec[[as.numeric(numberRe())]]
      tryCatch(load(paste0("www/rda/", df, ".rda")), error = function(e) {
        ""
      }, warning = function(w) {
        ""
      })
      values[["df1Mark"]] <- get(df)
      values[["df1MarkName"]] <- ifelse(df == "emptydata.frame", "no", df)
    }
  }
)

observeEvent(input$file1Mark, {
  inFile <- input[["file1Mark"]]

  if (grepl(".csv", inFile$datapath)) {
    values[["df1Mark"]] <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(".rds", inFile$datapath)) {
    values[["df1Mark"]] <- readRDS(inFile$datapath)
  }

  values[["df1MarkName"]] <- inFile$name
  values[["df1MarkOrigin"]] <- "file was loaded"
})

#
#   update mark style data.frame
#

observeEvent(numberRe(),
  ignoreInit = T,
  {
    validate(need(try(numberRe() > 0 & numberRe() <= maxEx), "not ready"))
    if (values[["canButton"]]) {
      values[["df1MStyle"]] <- data.frame()
      df <- values[["paramVec"]]$dfMarkColorVec[[as.numeric(numberRe())]]
      tryCatch(load(paste0("www/rda/", df, ".rda")), error = function(e) {
        ""
      }, warning = function(w) {
        ""
      })
      values[["df1MStyle"]] <- get(df)
      values[["df1MStyleName"]] <- ifelse(df == "emptydata.frame", "no", df)
    }
  }
)

observeEvent(input$file1MStyle, {
  inFile <- input[["file1MStyle"]]

  if (grepl(".csv", inFile$datapath)) {
    values[["df1MStyle"]] <- read.csv(inFile$datapath, header = TRUE, stringsAsFactors = FALSE)
  } else if (grepl(".rds", inFile$datapath)) {
    values[["df1MStyle"]] <- readRDS(inFile$datapath)
  }

  values[["df1MStyleName"]] <- inFile$name
  values[["df1MStyleOrigin"]] <- "file was loaded"
})

#
#   keep inputs notes
#

observeEvent(input$outNotes, {
  if (!is.null(input$outNotes)) {
    tryCatch(df <- hot_to_r(input$outNotes), error = function(e) {
      data.frame()
    })
    if (!is.function(df)) {
      values[["notes"]] <- df
    }
  }
})

#
#   keep inputs leftNotes
#

observeEvent(input$outleftNotes, {
  if (!is.null(input$outleftNotes)) {
    tryCatch(df <- hot_to_r(input$outleftNotes), error = function(e) {
      data.frame()
    })
    if (!is.function(df)) {
      values[["leftNotes"]] <- df
    }
  }
})

#
#   keep inputs leftNotesUp
#

observeEvent(input$outleftNotesUp, {
  if (!is.null(input$outleftNotesUp)) {
    tryCatch(df <- hot_to_r(input$outleftNotesUp), error = function(e) {
      data.frame()
    })
    if (!is.function(df)) {
      values[["leftNotesUp"]] <- df
    }
  }
})

#
#   keep inputs chr
#

observeEvent(input$out, {
  if (!is.null(input$out)) {
    tryCatch(df <- hot_to_r(input$out), error = function(e) {
      data.frame()
    })
    if (!is.function(df)) {
      values[["df1"]] <- df
    }
  }
})

#
#   keep inputs mark
#

observeEvent(input$outMark, {
  if (!is.null(input$outMark)) {
    tryCatch(df <- hot_to_r(input$outMark), error = function(e) {
      NA
      NA
    })
    if (!is.function(df)) {
      values[["df1Mark"]] <- df
    }
  }
})


#
#   keep inputs mark style
#

observeEvent(input$outMStyle, {
  tryCatch(df <- hot_to_r(input$outMStyle), error = function(e) {
    data.frame()
  })
  if (!is.function(df)) {
    values[["df1MStyle"]] <- df
  }
  # }
})

#
#   update filenames among analogous inputs
#

observeEvent(input$notesfilename, {
  updateTextInput(
    session = session,
    inputId = "notesfilename2",
    value = isolate(input$notesfilename)
  )
})

observeEvent(input[["notesfilename2"]], {
  updateTextInput(
    session = session,
    inputId = "notesfilename",
    value = isolate(input$notesfilename2)
  )
})

observeEvent(input$file1Notes, {
  updateTextInput(
    session = session,
    inputId = "notesfilename",
    value = isolate(values[["notesName"]])
  )
})


observeEvent(input$file1Notes, {
  updateTextInput(
    session = session,
    inputId = "notesfilename2",
    value = isolate(values[["notesName"]])
  )
})

#############

observeEvent(input$leftNotesfilename, {
  updateTextInput(
    session = session,
    inputId = "leftNotesfilename2",
    value = isolate(input$leftNotesfilename)
  )
})

observeEvent(input[["leftNotesfilename2"]], {
  updateTextInput(
    session = session,
    inputId = "leftNotesfilename",
    value = isolate(input$leftNotesfilename2)
  )
})

observeEvent(input$file1leftNotes, {
  updateTextInput(
    session = session,
    inputId = "leftNotesfilename",
    value = isolate(sub(".csv", "", values[["leftNotesName"]]))
  )
})

observeEvent(input$file1leftNotes, {
  updateTextInput(
    session = session,
    inputId = "leftNotesfilename2",
    value = isolate(sub(".csv", "", values[["leftNotesName"]]))
  )
})

###############

observeEvent(input$leftNotesUpfilename, {
  updateTextInput(
    session = session,
    inputId = "leftNotesUpfilename2",
    value = isolate(input$leftNotesUpfilename)
  )
})

observeEvent(input[["leftNotesUpfilename2"]], {
  updateTextInput(
    session = session,
    inputId = "leftNotesUpfilename",
    value = isolate(input$leftNotesUpfilename2)
  )
})

observeEvent(input$file1leftNotesUp, {
  updateTextInput(
    session = session,
    inputId = "leftNotesUpfilename",
    value = isolate(sub(".csv", "", values[["leftNotesUpName"]]))
  )
})

observeEvent(input$file1leftNotesUp, {
  updateTextInput(
    session = session,
    inputId = "leftNotesUpfilename2",
    value = isolate(sub(".csv", "", values[["leftNotesUpName"]]))
  )
})

#############
observeEvent(input$chrfilename, {
  updateTextInput(
    session = session,
    inputId = "chrfilename2",
    value = isolate(input$chrfilename)
  )
})

observeEvent(input[["chrfilename2"]], {
  updateTextInput(
    session = session,
    inputId = "chrfilename",
    value = isolate(input$chrfilename2)
  )
})

observeEvent(input$file1, {
  updateTextInput(
    session = session,
    inputId = "chrfilename",
    value = isolate(sub(".csv", "", values[["df1Name"]]))
  )
})

observeEvent(input$file1, {
  updateTextInput(
    session = session,
    inputId = "chrfilename2",
    value = isolate(sub(".csv", "", values[["df1Name"]]))
  )
})

################## 33

observeEvent(input$markfilename, {
  updateTextInput(
    session = session,
    inputId = "markfilename2",
    value = isolate(input$markfilename)
  )
})

observeEvent(input[["markfilename2"]], {
  updateTextInput(
    session = session,
    inputId = "markfilename",
    value = isolate(input$markfilename2)
  )
})

observeEvent(input[["file1Mark"]], {
  updateTextInput(
    session = session,
    inputId = "markfilename",
    value = isolate(sub(".csv", "", values[["df1MarkName"]]))
  )
})

observeEvent(input[["file1Mark"]], {
  updateTextInput(
    session = session,
    inputId = "markfilename2",
    value = isolate(sub(".csv", "", values[["df1MarkName"]]))
  )
})

##########

observeEvent(input$MStylefilename, {
  updateTextInput(
    session = session,
    inputId = "MStylefilename2",
    value = isolate(input$MStylefilename)
  )
})

observeEvent(input[["MStylefilename2"]], {
  updateTextInput(
    session = session,
    inputId = "MStylefilename",
    value = isolate(input$MStylefilename2)
  )
})

observeEvent(input[["file1MStyle"]], {
  updateTextInput(
    session = session,
    inputId = "MStylefilename",
    value = isolate(sub(".csv", "", values[["df1MStyle"]]))
  )
})

observeEvent(input[["file1MStyle"]], {
  updateTextInput(
    session = session,
    inputId = "MStylefilename2",
    value = isolate(sub(".csv", "", values[["df1MStyle"]]))
  )
})
