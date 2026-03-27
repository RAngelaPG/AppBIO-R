observeEvent(input$ARCIbutton, {
  values[["ARCImessage"]] <- ""
  if (!invalid(values[["df1"]])) {
    values[["ARCI"]] <- armRatioCI(values[["df1"]], rnumeric = TRUE)
    if (nrow(values[["ARCI"]]) == 0) {
      values[["ARCImessage"]] <- "Unable, check mandatory columns shortArmSize, longArmSize"
    }
  }
})

observeEvent(input$AA2button, {
  # values[["ARCImessage"]]<-""
  if (!invalid(values[["df1"]])) {
    res2 <- asymmetry(values[["df1"]], asDf = TRUE)
    values[["AA2"]] <- res2
  }
})

observeEvent(input$posCalcbutton, {
  # values[["ARCImessage"]]<-""
  if (!invalid(values[["df1"]]) & !invalid(values[["df1Mark"]])) {
    dfMark <- fillMarkInfo(values[["df1Mark"]], values[["df1"]], markDistType = input$markDistType, origin = input$origin)
    res2 <- posCalc(dfMark, values[["df1"]], result = "data.frame", bToRemove = values[["bToRemove"]])
    values[["posCalc"]] <- res2
  }
})

observeEvent(input$perMarkbutton, {
  # values[["ARCImessage"]]<-""
  if (!invalid(values[["df1"]]) & !invalid(values[["df1Mark"]])) {
    dfMark <- fillMarkInfo(values[["df1Mark"]], values[["df1"]], markDistType = input$markDistType, origin = input$origin)
    res2 <- perMark(dfMark, values[["df1"]], result = "data.frame", bToRemove = values[["bToRemove"]])
    values[["perMark"]] <- res2
  }
})
