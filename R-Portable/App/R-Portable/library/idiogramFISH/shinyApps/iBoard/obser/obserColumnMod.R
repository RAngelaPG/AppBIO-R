# removeFactorMark  update reactivs

observeEvent(input$removeFactorMark, {
  req(inherits(values[["df1Mark"]], "data.frame"))
  req(nrow(values[["df1Mark"]]) > 0)
  values[["df1Mark"]] <- tryCatch(idiogramFISH:::makeNumCols(values[["df1Mark"]]),
    error = function(e) {
      data.frame()
    }
  )
})

observeEvent(input$removeFactorChr, {
  req(inherits(values[["df1"]], "data.frame"))
  req(nrow(values[["df1"]]) > 0)
  values[["df1"]] <- tryCatch(idiogramFISH:::makeNumCols(values[["df1"]]),
    error = function(e) {
      data.frame()
    }
  )
})


#
#   add or remove columns chr
#

observeEvent(input$addColumn, {
  final <- values[["df1"]]
  if (inherits(final, "data.frame")) {
    if (nrow(final) > 0) {
      if (input$colType == "character") {
        final[, input$col] <- as.character(NA)
      } else {
        final[, input$col] <- as.numeric(NA)
      }
    } else {
      final <- data.frame(
        OTU = "1",
        chrName = "chrName",
        shortArmSize = 1,
        longArmSize = 1,
        chrSize = as.numeric(NA)
      )
    }
  } else {
    final <- data.frame(
      OTU = "1",
      chrName = "chrName",
      shortArmSize = 1,
      longArmSize = 1,
      chrSize = as.numeric(NA)
    )
  }
  values[["df1"]] <- final
})

observeEvent(input$removeColumn, {
  final <- values[["df1"]]
  final[, input$col] <- NULL
  values[["df1"]] <- final
})

#
#   notes add remove
#

observeEvent(input$addColumnNotes, {
  final <- values[["notes"]]

  if (inherits(final, "data.frame")) {
    if (nrow(final) > 0) {
      if (input$colTypeNotes == "character") {
        final[, input$colNotes] <- as.character(NA)
      } else {
        final[, input$colNotes] <- as.numeric(NA)
      }
    } else { # if no rows
      if (nrow(values[["df1"]]) > 0) {
        if (!"OTU" %in% colnames(values[["df1"]])) {
          values[["df1"]]$OTU <- as.character(1)
        }

        final <- data.frame(
          OTU = unique(values[["df1"]]$OTU),
          note = "note"
        )
      } else {
        final <- data.frame(
          OTU = "1",
          note = "note"
        )
      }
    }
  } else {
    if (nrow(values[["df1"]]) > 0) {
      if (!"OTU" %in% colnames(values[["df1"]])) {
        values[["df1"]]$OTU <- as.character(1)
      }

      final <- data.frame(
        OTU = unique(values[["df1"]]$OTU),
        note = "note"
      )
    } else {
      final <- data.frame(
        OTU = "1",
        note = "note"
      )
    }
  }
  values[["notes"]] <- final
})

observeEvent(input$removeColumnNotes, {
  final <- values[["notes"]]
  final[, input$colNotes] <- NULL
  values[["notes"]] <- final
})

####

#
#   leftNotes add remove
#

observeEvent(input$addColumnleftNotes, {
  final <- values[["leftNotes"]]

  if (inherits(final, "data.frame")) {
    if (nrow(final) > 0) {
      if (input$colTypeleftNotes == "character") {
        final[, input$colleftNotes] <- as.character(NA)
      } else {
        final[, input$colleftNotes] <- as.numeric(NA)
      }
    } else { # if no rows
      if (nrow(values[["df1"]]) > 0) {
        if (!"OTU" %in% colnames(values[["df1"]])) {
          values[["df1"]]$OTU <- as.character(1)
        }

        final <- data.frame(
          OTU = unique(values[["df1"]]$OTU),
          note = "note"
        )
      } else {
        final <- data.frame(
          OTU = "1",
          note = "note"
        )
      }
    }
  } else {
    if (nrow(values[["df1"]]) > 0) {
      if (!"OTU" %in% colnames(values[["df1"]])) {
        values[["df1"]]$OTU <- as.character(1)
      }

      final <- data.frame(
        OTU = unique(values[["df1"]]$OTU),
        note = "note"
      )
    } else {
      final <- data.frame(
        OTU = "1",
        note = "note"
      )
    }
  }
  values[["leftNotes"]] <- final
})

observeEvent(input$removeColumnleftNotes, {
  final <- values[["leftNotes"]]
  final[, input$colleftNotes] <- NULL
  values[["leftNotes"]] <- final
})

####
#
#   leftNotesUp add remove
#

observeEvent(input$addColumnleftNotesUp, {
  final <- values[["leftNotesUp"]]

  if (inherits(final, "data.frame")) {
    if (nrow(final) > 0) {
      if (input$colTypeleftNotesUp == "character") {
        final[, input$colleftNotesUp] <- as.character(NA)
      } else {
        final[, input$colleftNotesUp] <- as.numeric(NA)
      }
    } else { # if no rows
      if (nrow(values[["df1"]]) > 0) {
        if (!"OTU" %in% colnames(values[["df1"]])) {
          values[["df1"]]$OTU <- as.character(1)
        }

        final <- data.frame(
          OTU = unique(values[["df1"]]$OTU),
          note = "note"
        )
      } else {
        final <- data.frame(
          OTU = "1",
          note = "note"
        )
      }
    }
  } else {
    if (nrow(values[["df1"]]) > 0) {
      if (!"OTU" %in% colnames(values[["df1"]])) {
        values[["df1"]]$OTU <- as.character(1)
      }

      final <- data.frame(
        OTU = unique(values[["df1"]]$OTU),
        note = "note"
      )
    } else {
      final <- data.frame(
        OTU = "1",
        note = "note"
      )
    }
  }
  values[["leftNotesUp"]] <- final
})

observeEvent(input$removeColumnleftNotesUp, {
  final <- values[["leftNotesUp"]]
  final[, input$colleftNotesUp] <- NULL
  values[["leftNotesUp"]] <- final
})

####



#
#   marks add remove
#

observeEvent(input$addColumnMark, {
  final <- values[["df1Mark"]]
  if (inherits(final, "data.frame")) {
    if (nrow(final) > 0) {
      if (input$colTypeMark == "character") {
        final[, input$colMark] <- as.character(NA)
      } else {
        final[, input$colMark] <- as.numeric(NA)
      }
    } else { # if no rows

      final <- data.frame(
        OTU = "1",
        chrName = "chrName",
        markName = "markName",
        chrRegion = "p",
        markPos = 0.2,
        markSize = 0.2
      )

      if ("OTU" %in% colnames(values[["df1"]])) {
        final$OTU <- values[["df1"]]$OTU[1] # as.character(NA)
      }
    }
  } else {
    final <- data.frame(
      OTU = "1",
      chrName = "chrName",
      markName = "markName",
      chrRegion = "p",
      markPos = 0.2,
      markSize = 0.2
    )

    if ("OTU" %in% colnames(values[["df1"]])) {
      final$OTU <- values[["df1"]]$OTU[1] # final$OTU <-as.character(NA)
    }
  }
  values[["df1Mark"]] <- final
})


observeEvent(input$removeColumnMark, {
  final <- values[["df1Mark"]]
  final[, input$colMark] <- NULL
  values[["df1Mark"]] <- final
})

#
#   add or remove columns mark style
#

observeEvent(input$addColumnMStyle, {
  final <- values[["df1MStyle"]]
  if (inherits(final, "data.frame")) {
    if (nrow(final) > 0) {
      if (input$colTypeMStyle == "character") {
        final[, input$colMStyle] <- as.character(NA)
      } else {
        final[, input$colMStyle] <- as.numeric(NA)
      }
    } else {
      final <- data.frame(markName = "markName", markColor = "green", style = "dots")
    }
  } else {
    final <- data.frame(markName = "markName", markColor = "green", style = "dots")
  }
  values[["df1MStyle"]] <- final
})

observeEvent(input$removeColumnMStyle, {
  final <- values[["df1MStyle"]]
  final[, input$colMStyle] <- NULL
  values[["df1MStyle"]] <- final
})
