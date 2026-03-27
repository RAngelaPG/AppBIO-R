# outputNuccore
# searchUI  wellPanel

output$searchUI <- renderUI({
  conditionalPanel(
    condition = "input.fileInsteadNcbi == false",
    wellPanel(
      fluidRow(
        column(
          6,
          h4(tags$strong("String to search")),
          tags$span("with ", tags$code("rentrez::entrez_search")),
          h5("db = 'nuccore'")
        ),
        column(
          6,
          numericInput("maxNum",
            "max. number of results",
            value = 250,
            min = 1,
            max = 2500
          )
        )
      ),
      textInput("term",
        "",
        value = "",
        width = "100%"
      ),
      uiOutput("termButtonUI")
    )
  )
})


output$saveEntrezButton <- downloadHandler(
  filename = function() {
    values[["entrezFile"]]
  },
  content = function(file) {
    writeLines(values[["rentrezFetch"]], file)
  }
)


output$resultNumber <- renderText({
  validate(
    need(try(!is.na(values[["titles_number"]])), "")
  )
  paste(values[["titles_number"]], "results (max.", maxNumReac(), ")")
})

output$statsDF <- renderUI({
  validate(
    need(try({
      values[["annotated"]] == TRUE
    }), "")
  )
  validate(
    need(try(inherits(values[["geneMarkDFOrig"]], "data.frame")), "")
  )
  validate(
    need(try(inherits(values[["geneChrDF"]], "data.frame")), "")
  )
  wellPanel(
    h5(tags$strong("data.frame summary")),
    tagList(
      # tags$span(paste(nrow(values[["geneChrDF"]]),"rows in chr. d.f.") )
      # ,br()
      # ,
      tags$span(paste(nrow(values[["geneMarkDFOrig"]]), "rows in marks' d.f."))
    )
  )
})


output$chooseFileUI <- renderUI({
  conditionalPanel(
    condition = "input.fileInsteadNcbi == true",
    wellPanel(
      h4(tags$strong("Upload text file")),
      helpText("File downloaded with: rentrez::entrez_fetch(db='nuccore',
                                     id= [ID],
                                     rettype='gbwithparts',
                                     retmode = 'text')"),
      uiOutput("nucFileUI")
    ) # wP
  )
})

output$nucFileUI <- renderUI({
  input$modifyMarksButton
  input$nucPresetButton
  fileInput("nucFile", "3a. Choose:",
    accept = c(
      "text/plain",
      ".gb"
    )
  )
})

output$titleSelectUI <- renderUI({
  if (length(input$termButton == 0)) {
    if (input$termButton == 0) {
      return("Press 2. Search")
    } else {
      validate(
        need(try({
          values[["rentrezPkg"]] == TRUE
        }), values[["renMiss"]])
      )

      validate(
        need(try({
          values[["searchStatus"]]
        }), values[["errorMessage"]])
      )

      validate(
        need(try({
          !is.na(values[["entrez_titles"]])
          !is.na(values[["entrez_selected"]])
        }), "Press Search, and wait")
      )
      conditionalPanel(
        condition = "input.fileInsteadNcbi == false",
        wellPanel(
          h4(
            "Results of", tags$code("rentrez::entrez_search"), "and",
            tags$code("rentrez::entrez_summary")
          ),
          span(textOutput("resultNumber"), style = "font-size:15px"),
          selectizeInput("titleSelect",
            "Select match",
            choices = values[["named_entrez_titles"]],
            selected = values[["entrez_selected"]],
            options = list(maxOptions = 2500)
          ),
          actionButton("button3Download", tags$span(tags$strong("3b. Download with"), tags$code("rentrez::entrez_fetch"))),
          br(),
          uiOutput("saveEntrezUI")
        )
      )
    }
  } # termB not 0
})

output$saveEntrezUI <- renderUI({
  validate(
    need(
      try({
        class(values[["rentrezFetch"]]) == "character"
      }),
      "Press 3b."
    )
  )
  downloadButton("saveEntrezButton", tags$span(tags$strong("Save sequence (optional)")))
})


output$authorsUI <- renderUI({
  # if(length(input$termButton==0))
  #   if(input$termButton==0) {
  #     return("")
  #   } else {

  # validate(
  #   need(try({
  #     !is.na(values[["authors"]] )
  #   }
  #   ), "")
  # )

  wellPanel(
    h4(strong("Authors of sequence:")),
    htmlOutput("authors")
  )
  # }
})

output$termButtonUI <- renderUI({
  validate(
    need(
      try({
        input$term != ""
      }),
      "Fill manually or load preset"
    )
  )
  actionButton("termButton", tags$span(tags$strong("2. Search with"), tags$code("rentrez")))
})

output$loadUI <- renderUI({
  validate(need(try(values[["geneChrDF"]]), "Waiting for data.frames"))
  validate(need(try(inherits(values[["geneChrDF"]], "data.frame")), "still not ready chr. d.f.1"))
  validate(need(try(values[["geneMarkDF"]]), "Waiting for data.frames"))
  validate(need(try(inherits(values[["geneMarkDF"]], "data.frame")), "still not ready mark d.f.2"))
  validate(need(try(values[["markStyleDF"]]), "Waiting for data.frames"))
  validate(need(try(inherits(values[["markStyleDF"]], "data.frame")), "still not ready mark style d.f.3"))

  wellPanel(
    h4(strong("Data.frames for plotting")),
    h5(strong("Chr. marks and marks' style ")),
    actionButton("loadDFbutton", tags$strong("6. Load to data.frame page (and wait for plot)"))
  )
})

output$workflowUI <- renderUI({
  conditionalPanel(
    condition = "input.showWorkflow == true",
    wellPanel(
      h4(strong("Workflow")),
      tags$p("After each button click, a pop-up will indicate you have to wait;
              the pop-up will close automatically;
              then, the following button will be available"),
      h5(strong("Mandatory steps - buttons")),
      tags$ul( # tags$li("2: Search")
        tags$li("3a: Choose file or 2: Search and 3b. download data"),
        tags$li("4: make data.frames"),
        tags$li("6: load data.frames"),
        tags$li("Go to previous page for plot")
      ),
      br(),
      h5(strong("Additional waiting time")),
      tags$p("After button 6. pop-up, some additional time is required"),
      br(),
      tags$p("Wait some seconds for a plasmid plot (pop-up + plotting time)"),
      tags$p("For a chromosome, one minute (pop-ups) plus 30 seconds to plot"),
      br(),
      h5(strong("Where to go?")),
      tags$p("To check if data.frames have been loaded,
               go to the data.frames page; or wait in the plot page")
    ) # wP
  )
})

output$presetUI <- renderUI({
  wellPanel(
    radioButtons("nucPreset",
      "9. Circular plots: (nuc.)",
      inline = FALSE,
      choices = c(
        "9.3 plasmid preset" = 9,
        "9.4 bac. chrom. preset" = 10
      )
    )
    # c("Plasmid"=9,"Bac. Chromosome"=10) )
    , actionButton("nucPresetButton", tags$strong("1. Load preset"))
  )
})


output$markColumnUI <- renderUI({
  wellPanel(
    fluidRow(
      column(
        6,
        h4(tags$strong("Marks")),
        checkboxInput(
          "addSTARTPos", "add START",
          paramValues$addSTARTPosDefault
        ),
        div(
          title = "recommended because marks with same name have same style and color",
          checkboxInput(
            "makeUnique", "unique mark names",
            paramValues$makeUniqueDefault
          )
        ),
        radioButtons(
          "nucMarkStyle", "Style",
          c(
            "Arrows (upArrow/downArrow)" = "Arrows"
            # ,"Square (square/squareLeft)"="Square"
            , "cM (cm/cMLeft)" = "cM"
          ),
          paramValues$nucMarkStyleDefault
        ),
        h5(strong("Change orientation")),
        checkboxInput(
          "mirror", "As complement",
          paramValues$mirrorDefault
        ),
        h5(strong("Pseudogenes")),
        radioButtons("pseudo", "Show:", c(
          "only pseudo." = "onlyPseudo",
          "only not-pseudo." = "removePseudo",
          "all" = "all"
        ),
        selected = paramValues$pseudoDefault
        ),
        div(
          title = "replace locus tag",
          h5(strong("Replace names with: (when available)"))
        ),
        div(
          title = "use gene name instead of locus tag",
          checkboxInput("useGeneNames", "gene", paramValues$useGeneNamesDefault)
        ),
        div(
          title = "use regulatory class instead of locus tag",
          checkboxInput("useRCNames", "regulatory class", paramValues$useRCNamesDefault)
        )
      ),
      column(
        6,
        h4(strong("Separate mark names in columns")),
        splitLayout(
          div(
            title = "colNumber: number of columns",
            numericInput("colNumber",
              HTML(paste("Num. of", "columns", sep = "<br/>")),
              paramValues$colNumberDefault,
              min = 1,
              max = 10,
              step = 1
            )
          ),
          div(
            title = "protrudingInt: numeric, spacing of columns in terms of width of chr. percent 1 = 100%. Defaults to 0.5",
            numericInput("protrudingInt",
              HTML(paste("column", "spacing", sep = "<br/>")),
              paramValues$protrudingIntDefault,
              min = 0.1,
              max = 6,
              step = 0.1
            )
          ) # div
        ) # sL
        , div(
          title = "amountofSpaces, number of spaces for each column",
          numericInput("amountofSpaces",
            "spaces:",
            paramValues$amountofSpacesDefault,
            min = 1,
            max = 50,
            step = 1
          )
        ),
        div(
          title = 'markType, use c("downArrow","upArrow","cM","cMLeft") or a subset',
          checkboxGroupInput("markType",
            "types of mark to arrange in columns",
            c("downArrow", "upArrow", "cM", "cMLeft"),
            selected = paramValues$markTypeDefault
          )
        ),
        div(
          style = "margin-bottom:-5px",
          title = "`mycolors`: optional, character vector with colors' names, which are associated automatically with marks according to their order in the data.frame of position of marks. See this ordering with `unique(dfMarkPos$markName)`. Argument example: `mycolors = c(\"red\",\"chartreuse3\",\"dodgerblue\")`. Not mandatory for plotting marks, package has default colors.
          ",
          textInput(
            "mycolors2",
            HTML(paste(
              "Marks' colors",
              tags$p("comma separated", style = "font-size: 80%;")
            )),
            paramValues$mycolors2Default
          )
        ),
        div(
          title = "when selected colors will differentiate features not marks",
          checkboxInput(
            "colorFeature", "color for features",
            paramValues$colorFeatureDefault
          )
        )
      )
    ),
    uiOutput("namesColumnButtonUI")
  )
})

output$namesColumnButtonUI <- renderUI({
  validate(
    need(try({
      input$makeDFsButton
    }), "IMPORTANT! Press (4. Make data.frames)")
  )

  validate(need(try(values[["geneChrDF"]]), "Waiting for data.frames"))
  validate(need(try(inherits(values[["geneChrDF"]], "data.frame")), "still not ready chr. d.f.1"))
  validate(need(try(values[["geneMarkDFOrig"]]), "Waiting for data.frames"))
  validate(need(try(inherits(values[["geneMarkDFOrig"]], "data.frame")), "still not ready mark d.f.2"))

  if (input$makeDFsButton == 0) {
    return("IMPORTANT! Press (4. Make data.frames)")
  }
  actionButton("modifyMarksButton", tags$strong("5. Modify marks"))
})


output$fetchSelectButtonUI <- renderUI({
  validate(
    need(try({
      !is.na(values[["button3ab"]])
    }), "IMPORTANT! load (3.) and find first")
  )
  actionButton("makeDFsButton", tags$strong("4. Make data.frames"))
})

output$authors <- renderUI({
  return(values[["authors"]])
})

output$fetchSelectUI <- renderUI({
  if (values[["termButtonVal"]] == 0) {
    return("")
  } else {
    validate(
      need(try({
        values[["annotated"]] == TRUE
      }), "Sorry, not annotated genome, nothing to do")
    )
    validate(
      need(try({
        values[["searchStatus"]] == TRUE
      }), "")
    )

    validate(
      need(try({
        !is.na(values[["button3ab"]])
      }), "IMPORTANT! Press 3. ")
    )
    wellPanel(
      h4(tags$strong("Features"), "gotten with", tags$code("rentrez::entrez_fetch")),
      tags$div(
        align = "left",
        class = "multicol",
        checkboxGroupInput("fetchSelect",
          "Select type of sequence",
          choiceNames  =  values[["names_fetch_list"]],
          choiceValues =  1:length(values[["names_fetch_list"]]),
          selected = values[["fetch_list_selected"]]
        )
      ),
      uiOutput("fetchSelectButtonUI")
      # ,actionButton("makeDFsButton",tags$strong("4. Make data.frames") )
    ) # wP
  }
})
