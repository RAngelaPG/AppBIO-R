observe({
  swapI <- tryCatch(unlist(strsplit(input$chrNamesToSwap, ",")),
    error = function(e) {
      NULL
    }
  )

  values[["chrNamesToSwap"]] <- swapI
})

#
#   to reactive
#

# observeEvent(input$exampleButton,{
#   values[["dfList"]] <- list()
# })

observeEvent(
  {
    input$swapButton
  },
  {
    chrNamesToSwapPar <- if (length(values[["chrNamesToSwap"]]) == 0) {
      ""
    } else {
      values[["chrNamesToSwap"]]
    }

    tryCatch(dfList <- swapChrRegionDfSizeAndMarks(values[["df1"]], values[["df1Mark"]], chrNamesToSwapPar),
      error = function(e) {
        "waiting for d.f."
      }
    )

    tryCatch(values[["dfList"]] <- dfList,
      error = function(e) {
        "waiting for d.f."
      }
    )

    if (inherits(values[["dfList"]]$dfChrSize, "data.frame")) {
      df1 <- values[["dfList"]]$dfChrSize
      values[["df1"]] <- df1
      df1Mark <- values[["dfList"]]$dfMarkPos
      values[["df1Mark"]] <- df1Mark
    }
  }
)
