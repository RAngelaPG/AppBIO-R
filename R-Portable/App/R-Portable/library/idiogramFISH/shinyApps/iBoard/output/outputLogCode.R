# output log code

output$log <- renderText({
  rawText <- filenameR()
  validate(need(try(readLines(rawText), TRUE), message = FALSE))

  replacedText <- paste(tryCatch(readLines(rawText),
    error = function(e) {
      "nothing to write"
    },
    warning = function(w) {
      "nothing to write"
    }
  ),
  collapse = "\n"
  )
  values[["log"]] <- replacedText
})

output$code <- renderText({
  return(values[["strFun"]])
})

output$logpanel <- renderUI({
  fluidRow(
    column(
      width = 6,
      wellPanel(
        uiOutput("cliplog"),
        h2("log"),
        tags$style(type = "text/css", "#log {background-color: #FFFFFF; color: black;}"),
        verbatimTextOutput("log")
      )
    )
  )
})

output$clip <- renderUI({
  rclipButton("clipbtn", "Copy code", values[["strFun"]], icon = icon("clipboard"))
})

output$cliplog <- renderUI({
  rclipButton("clipbtnlog", "Copy log", values[["log"]], icon = icon("clipboard"))
})
