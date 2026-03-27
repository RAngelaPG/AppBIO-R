output$currentExample <- renderText({
  paste(
    "Last example loaded:",
    names(exampleVec[grep(paste0("\\<", exampleVec[exampleVec == values[["number"]]], "\\>"), exampleVec)])
  )
})

output$wPUpPreset <- renderUI({
  validate(need(try(values[["presetStatus"]] == "Ready"), ""))

  if (current_lenR() > iniLen) {
    cn <- paste("Custom Preset", 1:(current_lenR() - iniLen))
    cv <- (iniLen + 1):current_lenR()

    wellPanel(
      h4("Custom Presets"),
      helpText("this will overwrite any data you have entered in other tabs"),
      radioButtons("upPreset",
        label = "",
        choices = setNames(cv, cn)
      ),
      actionButton("upPresetButton",
        tags$span(tags$strong("Load Presets")),
        class = "btn-success"
        # ,icon =
      )
    )
  }
})

output$examplepanel <- renderUI({
  fluidRow(
    column(
      width = 4,
      br(),
      wellPanel(
        h4("Choose an example"),
        helpText("numbered according to vignettes"),
        helpText("this will overwrite any data you have entered in other tabs"),
        radioGroupContainer(
          "exampleId",
          fluidRow(
            column(
              12,
              radioSubgroup("exampleId"
                # ,id="minimal_examples"
                ,
                label = "3. Minimal Examples:",
                choiceNames = c(
                  "3.1 monocentric k.",
                  "3.2 holocentric k.",
                  "3.3 mono. and holo.",
                  "3.4 circular m. & h. "
                ),
                choiceValues = c(1:2, 11:12),
                first = TRUE
              )
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "plotting_chr"
                ,
                label = "4. Plotting chr.:",
                choiceNames = c(
                  "4.2 monocen.",
                  "4.3 holocen.",
                  "4.3b holocentric k. mycolors",
                  "4.5 mono & holo in kar."
                ),
                choiceValues = c(13, 14, 5, 15)
              )
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "multiple_otus"
                ,
                label = "5. Multiple OTUs:",
                choiceNames = c(
                  "5.1 multiple mono. k.",
                  "5.2 multiple holo. k."
                ),
                choiceValues = 3:4
              )
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "changing_units"
                ,
                label = "6. Changing Units:",
                choiceNames = c(
                  "6.1 Using bases instead of micrometers",
                  "6.2 Using threshold to fix scale",
                  "6.3 Plot data in micrometers and bases",
                  "6.4 Use cM as units"
                ),
                choiceValues = 16:19
              )
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "gish"
                ,
                label = "7. GISH:",
                choiceNames = c(
                  "7.1 GISH",
                  "7.2 GISH holoc.",
                  "7.3a GISH citrus",
                  "7.3b GISH citrus"
                ),
                choiceValues = c(6, 7, 20, 21)
              )
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "gish"
                ,
                label = "8. Groups:",
                choiceNames = c(
                  "8.1 monoc.",
                  "8.2 holoc."
                ),
                choiceValues = c(22, 23)
              )
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "circular"
                ,
                label = "9. Circular plots:",
                choiceNames = c(
                  "9.1 Plot 1",
                  "9.2 Plot 2",
                  "9.3 plasmid preset",
                  "9.4 bac. chrom. preset"
                ),
                choiceValues = c(24, 8, 9, 10)
              ),
              helpText("9.3 & 9.4 redirect to Nucleotides")
            ),
            column(
              12,
              radioSubgroup("exampleId"
                # , "phylog"
                ,
                label = "10. For phylog:",
                choiceNames = c(
                  "10.2 Mono",
                  "10.3 Holo",
                  "10.4 mono & holo"
                ),
                choiceValues = c(25, 26, 27)
              ),
              helpText("See vignettes for the phylogeny")
            ),
            column(
              12,
              radioSubgroup("exampleId",
                label = tagList(
                  tags$span("11."),
                  tags$i("Citrus")
                ),
                choiceNames = c(
                  "11.1 C. maxima",
                  "11.2 C. reticulata",
                  "11.7 C. limon origin"
                ),
                choiceValues = c(28, 29, 30)
              )
            )
          )
        ),
        actionButton("exampleButton",
          tags$span(tags$strong("Load Example presets")),
          class = "btn-info"
        ),
        helpText("Click button to start")
      )
    ),
    column(
      width = 2,
      br(),
      wellPanel(
        h4("Upload Preset"),
        helpText("Presets can be downloaded in the 'Parameters' page"),
        fileInput("upFilePreset",
          "Choose .rds File",
          accept = c(
            ".rds"
          ),
          buttonLabel = list(icon("folder"))
        )
      ),
      uiOutput("wPUpPreset")
    ),
    column(
      width = 3,
      br(),
      wellPanel(
        h4("Jupyter notebooks"),
        helpText("Based on vignettes. Work online in colab or locally after downloading"),
        tags$table(
          tags$tr(
            tags$th(img(src = "colab-badge.svg")),
            tags$th(img(src = "GitHub-Mark-120px-plus.png", width = 18), HTML("GitHub&emsp;&emsp;")),
            tags$th("")
          ),
          add_row("03-minimal.ipynb", "3 Minimal examples"),
          add_row("04-plotting.ipynb", "4 Plotting chromosomes"),
          add_row("05-multiple.ipynb", "5 Multiple OTUs"),
          add_row("06-units.ipynb", "6 Changing units"),
          add_row("07-gish.ipynb", "7 GISH"),
          add_row("08-groups.ipynb", "8 Groups"),
          add_row("09-circular.ipynb", "9 Circular Plots"),
          add_row("10-phylogeny.ipynb", "10 Plotting alongside phylogeny"),
          add_row("11-citrushelp.ipynb", "11 Citrus"),
          add_row("12-human.ipynb", "12 Human Karyotype")
        )
      ),
      wellPanel(
        h4("Help"),
        helpText("Links to online vignettes and manual in the", tags$strong("About"), "page")
      )
    )
  )
})

output$ARCImessageOut <- renderText({
  return(values[["ARCImessage"]])
})

output$ARCImessageUI <- renderUI({
  req(values[["ARCImessage"]] != "")
  tagList(
    br(),
    verbatimTextOutput("ARCImessageOut")
  )
})

output$indicespanel <- renderUI({
  fluidRow(
    column(
      width = 2,
      br(),
      wellPanel(
        h5("Use function armRatioCI on"),
        h5("main data.frame (see data.frames page)"),
        helpText("reads a data.frame and produces AR (arm ratio), CI (centromeric index)
                       , Guerra (1986) and Levan et al. (1964) classifications."),
        actionButton("ARCIbutton",
          tags$span(tags$strong("armRatioCI")),
          class = "btn-info"
        ),
        uiOutput("ARCImessageUI"),
        uiOutput("downloadARCIUI"),
        br(),
        uiOutput("clipARCI")
      ),
      br(),
      wellPanel(
        h5("Use function asymmetry on"),
        h5("main data.frame"),
        helpText("calculates karyotype asymmetry A and A2."),
        actionButton("AA2button",
          tags$span(tags$strong("asymmetry")),
          class = "btn-info"
        ),
        uiOutput("downloadAA2UI"),
        br(),
        uiOutput("clipAA2")
      )
    ),
    column(
      7,
      br(),
      uiOutput("outARCIUI")
    ),
    column(
      2,
      br(),
      uiOutput("outAA2UI")
    )
  )
})

output$downloadARCIUI <- renderUI({
  req(nrow(values[["ARCI"]]) > 0)
  tagList(
    br(),
    h5("Save AR and CI"),
    textInput("ARCIfilename",
      "File name",
      value = "ARCIData",
      width = "50%"
    ),
    uiOutput("ARCIbuttontable")
  )
})

output$downloadAA2UI <- renderUI({
  req(nrow(values[["AA2"]]) > 0)
  tagList(
    br(),
    h5("Save A and A2"),
    textInput("AA2filename",
      "File name",
      value = "AA2Data",
      width = "50%"
    ),
    uiOutput("AA2buttontable")
  )
})

output$markspanel <- renderUI({
  fluidRow(
    column(
      2,
      br(),
      wellPanel(
        h5("Use function posCalc on"),
        h5("main data.frame (see data.frames page)"),
        helpText("calculates position of marks in fraction of (%) chromosome units (0-1)"),
        actionButton("posCalcbutton",
          tags$span(tags$strong("posCalc")),
          class = "btn-info"
        )
        # ,uiOutput("ARCImessageUI")
        , uiOutput("downloadposCalcUI"),
        br(),
        uiOutput("clipposCalc")
      ),
      br(),
      wellPanel(
        h5("Use function perMark on"),
        h5("main data.frame"),
        helpText("calculates fraction (%) of chromosome for each name of mark"),
        helpText("marks with location/chrRegion 'cen' don't have mark size"),
        actionButton("perMarkbutton",
          tags$span(tags$strong("perMark")),
          class = "btn-info"
        ),
        uiOutput("downloadperMarkUI"),
        br(),
        uiOutput("clipperMark")
      )
    ),
    column(
      5,
      br(),
      uiOutput("outposCalcUI")
    ),
    column(
      5,
      br(),
      uiOutput("outperMarkUI")
    )
  )
})

output$downloadposCalcUI <- renderUI({
  req(nrow(values[["posCalc"]]) > 0)
  tagList(
    br(),
    h5("Save posCalc"),
    textInput("posCalcfilename",
      "File name",
      value = "posCalcData",
      width = "50%"
    ),
    uiOutput("posCalcbuttontable")
  )
})

output$downloadperMarkUI <- renderUI({
  req(nrow(values[["perMark"]]) > 0)
  tagList(
    br(),
    h5("Save perMark"),
    textInput("perMarkfilename",
      "File name",
      value = "perMarkData",
      width = "50%"
    ),
    uiOutput("perMarkbuttontable")
  )
})


output$dfchrpanel <- renderUI({
  fluidRow(
    column(
      width = 2,
      br(),
      uiOutput("saveChrSize1UI"),
      wellPanel(
        h4("Upload data"),
        helpText("optional"),
        fileInput("file1", "Choose csv/rds File",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".rds"
          )
        )
      )
    ),
    column(
      width = 2,
      br(),
      wellPanel(
        h4("Add/Remove Row"),
        helpText(
          "",
          p("Right-click inside the table to delete/insert rows.")
        )
      ),
      wellPanel( # style="margin-top:-15px;padding: 5px 5px 5px 5px"
        h4("Columns"),
        helpText("Use unique names"),
        textInput("col",
          "Col. name",
          value = "newCol",
          width = "50%"
        ),
        radioButtons("colType", "Col. data type", c("character", "numeric")),
        actionButton(
          "addColumn",
          "Add"
        ),
        actionButton(
          "removeColumn",
          "Remove"
        ),
        br(),
        div(
          style = "margin-top:10px",
          actionButton("removeFactorChr", "Factors (dropdown) to character")
        )
      )
    ),
    column(
      8,
      br(),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            h4("swap arms"),
            textInput("chrNamesToSwap", "chr. names", "1,2"),
            actionButton("swapButton", "Swap!", icon = icon("sync"))
          )
        )
      ),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            h2("Chr. data data.frame"),
            helpText(paste(
              values[["df1Name"]],
              values[["df1Origin"]]
            ))
          )
        )
      ),
      fluidRow(
        column(
          11,
          uiOutput("outUI")
        )
      )
    )
  )
})

output$dfmarkpanel <- renderUI({
  fluidRow(
    column(
      width = 2,
      br(),
      uiOutput("saveMarkPos1UI"),
      wellPanel(
        h4("Upload data"),
        helpText("optional"),
        fileInput("file1Mark",
          "Choose csv/rds File",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv",
            ".rds"
          )
        )
      )
    ),
    column(
      width = 2,
      br(),
      wellPanel(
        h4("Add/Remove Row"),
        helpText(
          "",
          p("Right-click inside the table to delete/insert rows.")
        )
      ),
      wellPanel(
        h4("Columns"),
        helpText("Use unique names"),
        textInput("colMark",
          "Col. name",
          value = "newCol",
          width = "50%"
        ),
        radioButtons("colTypeMark", "Col. data type", c("character", "numeric")),
        actionButton(
          "addColumnMark",
          "Add"
        ),
        actionButton(
          "removeColumnMark",
          "Remove"
        ),
        br(),
        div(
          style = "margin-top:10px",
          actionButton("removeFactorMark", "Factors (dropdown) to character")
        )
      )
    ),
    column(
      8,
      br(),
      fluidRow(
        column(
          width = 6,
          wellPanel(
            h2("Mark pos. data.frame"),
            helpText(paste(values[["df1MarkName"]], values[["df1MarkOrigin"]]))
          ),
        )
      ),
      fluidRow(
        column(
          11,
          uiOutput("outMarkUI")
        )
      )
    )
  )
})

output$dfMStylepanel <- renderUI({
  fluidRow(
    column(
      width = 2,
      br(),
      uiOutput("saveMarkStyle1UI"),
      wellPanel(
        h4("Upload data"),
        helpText("optional"),
        fileInput("file1MStyle",
          "Choose csv/rds File",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv", ".rds"
          )
        )
      )
    ),
    column(
      width = 2,
      br(),
      wellPanel(
        h4("Add/Remove Row"),
        helpText(
          "",
          p("Right-click inside the table to delete/insert rows.")
        )
      ),
      wellPanel(
        h4("Columns"),
        helpText("Use unique names"),
        textInput("colMStyle",
          "Col. name",
          value = "newCol",
          width = "50%"
        ),
        radioButtons("colTypeMStyle", "Col. data type", c("character", "numeric")),
        actionButton(
          "addColumnMStyle",
          "Add"
        ),
        actionButton(
          "removeColumnMStyle",
          "Remove"
        )
      )
    ),
    column(
      width = 5,
      br(),
      wellPanel(
        h2("Mark style data.frame"),
        helpText(paste(values[["df1MStyleName"]], values[["df1MStyleOrigin"]]))
      ),
      uiOutput("outMStyleUI")
    )
  )
})





output$dfnotespanel <- renderUI({
  tagList(
    fluidRow(
      column(
        width = 2,
        br(),
        uiOutput("saveNotes1UI"),
        wellPanel(
          h4("Upload data"),
          helpText("optional"),
          fileInput("file1Notes",
            "Choose csv/rds File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv", ".rds"
            )
          )
        )
      ),
      column(
        width = 2,
        br(),
        wellPanel(
          h4("Add/Remove Row"),
          helpText(
            "",
            p("Right-click inside the table to delete/insert rows.")
          )
        ),
        wellPanel(
          style = "margin-top:-15px;padding: 5px 5px 5px 5px",
          h4("Columns"),
          helpText("Use unique names"),
          splitLayout(
            textInput("colNotes",
              "Col. name",
              value = "newCol",
              width = "95%"
            ),
            radioButtons("colTypeNotes", "Col. data type", c("character", "numeric"))
          ),
          actionButton(
            "addColumnNotes",
            "Add"
          ),
          actionButton(
            "removeColumnNotes",
            "Remove"
          ),
          br()
        )
      ),
      column(
        8,
        br(),
        fluidRow(
          column(
            width = 6,
            wellPanel(
              h2("4. notes data.frame"),
              helpText(paste(values[["notesName"]], values[["notesOrigin"]])),
              helpText("text to be shown to the right of the karyotype")
            )
          )
        ),
        fluidRow(
          column(
            11,
            uiOutput("outnotesUI")
          )
        )
      )
    ),
    fluidRow(
      column(
        width = 2,
        br(),
        uiOutput("saveleftNotes1UI"),
        wellPanel(
          h4("Upload data"),
          helpText("optional"),
          fileInput("file1leftNotes",
            "Choose csv/rds File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv", ".rds"
            )
          )
        )
      ),
      column(
        width = 2,
        br(),
        wellPanel(
          h4("Add/Remove Row"),
          helpText(
            "",
            p("Right-click inside the table to delete/insert rows.")
          )
        ),
        wellPanel(
          style = "margin-top:-15px;padding: 5px 5px 5px 5px",
          h4("Columns"),
          helpText("Use unique names"),
          splitLayout(
            textInput("colleftNotes",
              "Col. name",
              value = "newCol",
              width = "95%"
            ),
            radioButtons("colTypeleftNotes", "Col. data type", c("character", "numeric"))
          ),
          actionButton(
            "addColumnleftNotes",
            "Add"
          ),
          actionButton(
            "removeColumnleftNotes",
            "Remove"
          ),
          br()
        )
      ),
      column(
        8,
        br(),
        fluidRow(
          column(
            width = 6,
            wellPanel(
              h2("5. leftNotes data.frame"),
              helpText(paste(values[["leftNotesName"]], values[["leftNotesOrigin"]])),
              helpText("text to be shown to the left of the karyotype")
            )
          )
        ),
        fluidRow(
          column(
            11,
            uiOutput("outleftNotesUI")
          )
        )
      )
    ),
    fluidRow(
      column(
        width = 2,
        br(),
        uiOutput("saveleftNotesUp1UI"),
        wellPanel(
          h4("Upload data"),
          helpText("optional"),
          fileInput("file1leftNotesUp",
            "Choose csv/rds File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv", ".rds"
            )
          )
        )
      ),
      column(
        width = 2,
        br(),
        wellPanel(
          h4("Add/Remove Row"),
          helpText(
            "",
            p("Right-click inside the table to delete/insert rows.")
          )
        ),
        wellPanel(
          style = "margin-top:-15px;padding: 5px 5px 5px 5px",
          h4("Columns"),
          helpText("Use unique names"),
          splitLayout(
            textInput("colleftNotesUp",
              "Col. name",
              value = "newCol",
              width = "95%"
            ),
            radioButtons("colTypeleftNotesUp", "Col. data type", c("character", "numeric"))
          ),
          actionButton(
            "addColumnleftNotesUp",
            "Add"
          ),
          actionButton(
            "removeColumnleftNotesUp",
            "Remove"
          ),
          br()
        )
      ),
      column(
        8,
        br(),
        fluidRow(
          column(
            width = 6,
            wellPanel(
              h2("6. leftNotesUp data.frame"),
              helpText(paste(values[["leftNotesUpName"]], values[["leftNotesUpOrigin"]])),
              helpText("text to be shown to the up-left of the karyotype")
            )
          )
        ),
        fluidRow(
          column(
            11,
            uiOutput("outleftNotesUpUI")
          )
        )
      )
    )
  )
})

output$circParam <- renderUI({
  wellPanel(
    h4("circular plot param."),
    splitLayout(
      div(
        title = chrLabelSpacingDesc,
        numericInput("chrLabelSpacing",
          "chr. name dist.",
          paramValues$chrLabelSpacingDefault,
          min = 0.25, max = 10, step = 0.25
        )
      ),
      div(
        title = labelSpacingDesc,
        numericInput("labelSpacing",
          "label spacing",
          paramValues$labelSpacingDefault,
          min = 0.25, max = 10, step = 0.05
        )
      )
    ),
    splitLayout(
      div(
        title = rotationDesc,
        numericInput("rotation",
          "Rotation",
          paramValues$rotationDefault,
          min = 0, max = 2 * pi, step = 0.05
        )
      ),
      div(
        style = "margin-top: 0px; margin-bottom:0px",
        title = labelOutwardsDesc,
        checkboxInput(
          "labelOutwards",
          HTML(paste("outwards", "projected", "label marks", sep = "<br/>")),
          paramValues$labelOutwardsDefault
        )
      )
    ),
    splitLayout(
      div(
        title = shrinkFactorDesc,
        numericInput("shrinkFactor",
          "shrink kar.",
          paramValues$shrinkFactorDefault,
          min = 0.2, max = 1, step = 0.05
        )
      ),
      div(
        title = radiusDesc,
        numericInput("radius",
          "radius",
          paramValues$radiusDefault,
          min = 0.25, max = 30, step = 0.05
        )
      )
    ),
    splitLayout(
      div(
        title = separFactorDesc,
        numericInput("separFactor",
          "separ kar.",
          paramValues$separFactorDefault,
          min = 0.5, max = 10, step = 0.5
        )
      ),
      div(
        title = circleCenterDesc,
        numericInput("circleCenter",
          "pos. center X-axis",
          paramValues$circleCenterDefault,
          min = -5, max = 15, step = 0.5
        )
      )
    )
  )
})

output$circParamOTU <- renderUI({
  wellPanel(
    h3("OTU name"),
    h5("(circ. plot)"),
    fluidRow(
      column(
        6,
        div(
          title = OTUplacingDesc,
          radioButtons("OTUplacing", "position",
            c("first", "number", "simple"),
            selected = paramValues$OTUplacingDefault
          )
        )
      ), column(
        6,
        div(
          title = OTUsrtDesc,
          numericInput("OTUsrt",
            "Angle",
            paramValues$OTUsrtDefault,
            min = -360, max = 360, step = 1
          )
        ),
        div(
          title = OTUjustifDesc,
          numericInput("OTUjustif",
            "Justif.",
            paramValues$OTUjustifDefault,
            min = 0, max = 1, step = 0.5
          )
        )
      )
    ),
    splitLayout(
      div(
        title = OTULabelSpacerxDesc,
        numericInput("OTULabelSpacerx",
          HTML(paste("position X-axis")),
          paramValues$OTULabelSpacerxDefault,
          min = -10, max = 10, step = .5
        )
      ),
      div(
        title = OTUlegendHeightDesc,
        numericInput("OTUlegendHeight",
          HTML(paste("vert. size")),
          paramValues$OTUlegendHeightDefault,
          min = 0.5, max = 10, step = 0.5
        )
      )
    )
  )
})

output$imageWidth <- renderUI({
  px <- session$clientData$output_idiogramPlot_width * (as.numeric(input$widFactor) / 100)
  px_formatted <- format(round(px, 2), nsmall = 2)
  px_formatted_80 <- format(round(px / 80, 2), nsmall = 2)
  HTML(paste(paste(px_formatted, "px"), paste(px_formatted_80, "in"), sep = "<br/>"))
})

output$imageHeight <- renderUI({
  px <- session$clientData$output_idiogramPlot_width * (as.numeric(input$heiFactor) * (as.numeric(input$widFactor) / 100))
  px_formatted <- format(round(px, 2), nsmall = 2)
  px_formatted_80 <- format(round(px / 80, 2), nsmall = 2)
  HTML(paste(paste(px_formatted, "px"), paste(px_formatted_80, "in"), sep = "<br/>"))
})

output$strpanel <- renderUI({
  fluidRow(
    column(
      width = 6,
      br(),
      wellPanel( # style = "padding: 0px 5px 0px 5px"
        splitLayout(
          cellWidths = c("33%", "66%"),
          div(
            style = "max-width:100px;",
            uiOutput("buttonPresets")
          ),
          helpText(presetUseText)
        )
      ),
      wellPanel(
        splitLayout(
          uiOutput("clip"),
          checkboxInput("keepDesc", "keep description", TRUE),
          checkboxInput("asFile", "plot to file (keeps size)", TRUE)
        ),
        splitLayout(
          div(
            style = "max-width:100px;",
            uiOutput("buttonScript")
          ),
          checkboxInput("keepDefault", "keep default values", FALSE),
          radioButtons("saveType", "download format", c("csv", "rds"), "csv", inline = T)
        ),
        tags$style(type = "text/css", "#code {background-color: #FFFFFF; color: black;}"),
        verbatimTextOutput("code")
      )
    ),
    column(
      width = 2,
      br(),
      uiOutput("saveChrSizeUI"),
      uiOutput("saveMarkPosUI"),
      uiOutput("saveMarkStyleUI")
    ),
    column(
      2,
      br(),
      uiOutput("saveNotesUI"),
      uiOutput("saveleftNotesUI"),
      uiOutput("saveleftNotesUpUI")
    )
  )
})

output$saveChrSize1UI <- renderUI({
  if (!invalid(values[["df1"]])) {
    wellPanel(
      h4("Save data"),
      helpText("optional"),
      textInput("chrfilename",
        "File name",
        value = "chrData",
        width = "50%"
      ),
      uiOutput("buttontable")
    )
  }
})

output$saveChrSizeUI <- renderUI({
  if (!invalid(values[["df1"]])) {
    wellPanel(
      h4("Save chr. data"),
      helpText(""),
      textInput("chrfilename2",
        "File name",
        value = "chrData",
        width = "50%"
      ),
      uiOutput("buttontable2")
    )
  }
})

output$saveMarkPos1UI <- renderUI({
  if (!invalid(values[["df1Mark"]])) {
    wellPanel(
      h4("Save data"),
      helpText("optional"),
      textInput("markfilename",
        "File name",
        value = "markData",
        width = "50%"
      ),
      uiOutput("buttontableMark")
    )
  }
})

output$saveMarkPosUI <- renderUI({
  if (!invalid(values[["df1Mark"]])) {
    wellPanel(
      h4("Save mark pos. data"),
      helpText(""),
      textInput("markfilename2",
        "File name",
        value = "markData",
        width = "50%"
      ),
      uiOutput("buttontableMark2")
    )
  }
})

output$saveMarkStyle1UI <- renderUI({
  if (!invalid(values[["df1MStyle"]])) {
    wellPanel(
      h4("Save data"),
      helpText("optional"),
      textInput("MStylefilename",
        "File name",
        value = "MStyleData",
        width = "50%"
      ),
      uiOutput("buttontableMStyle")
    )
  }
})

output$saveMarkStyleUI <- renderUI({
  if (!invalid(values[["df1MStyle"]])) {
    wellPanel(
      h4("Save mark style data"),
      helpText(""),
      textInput("MStylefilename2",
        "File name",
        value = "MStyleData",
        width = "50%"
      ),
      uiOutput("buttontableMStyle2")
    )
  }
})

output$saveNotes1UI <- renderUI({
  if (!invalid(values[["notes"]])) {
    wellPanel(
      style = "margin-bottom:0px;padding: 5px 5px 5px 5px",
      splitLayout(
        tagList(
          h4("Save data"),
          helpText("optional")
        ),
        textInput("notesfilename",
          "File name",
          value = "notes",
          width = "90%"
        )
      ),
      uiOutput("buttontablenotes")
    )
  }
})

output$saveNotesUI <- renderUI({
  if (!invalid(values[["notes"]])) {
    wellPanel(
      h4("Save notes data.frame"),
      helpText(""),
      textInput("notesfilename2",
        "File name",
        value = "notes",
        width = "50%"
      ),
      uiOutput("buttontablenotes2")
    )
  }
})


output$saveleftNotes1UI <- renderUI({
  if (!invalid(values[["leftNotes"]])) {
    wellPanel(
      style = "margin-bottom:0px;padding: 5px 5px 5px 5px",
      splitLayout(
        tagList(
          h4("Save data"),
          helpText("optional")
        ),
        textInput("leftNotesfilename",
          "File name",
          value = "leftNotes",
          width = "90%"
        )
      ),
      uiOutput("buttontableleftNotes")
    )
  }
})

output$saveleftNotesUI <- renderUI({
  if (!invalid(values[["leftNotes"]])) {
    wellPanel(
      h4("Save leftNotes data.frame"),
      helpText(""),
      textInput("leftNotesfilename2",
        "File name",
        value = "leftNotes",
        width = "50%"
      ),
      uiOutput("buttontableleftNotes2")
    )
  }
})

output$saveleftNotesUp1UI <- renderUI({
  if (!invalid(values[["leftNotesUp"]])) {
    wellPanel(
      style = "margin-bottom:0px;padding: 5px 5px 5px 5px",
      splitLayout(
        tagList(
          h4("Save data"),
          helpText("optional")
        ),
        textInput("leftNotesUpfilename",
          "File name",
          value = "leftNotesUp",
          width = "90%"
        )
      ),
      uiOutput("buttontableleftNotesUp")
    )
  }
})

output$saveleftNotesUpUI <- renderUI({
  if (!invalid(values[["leftNotesUp"]])) {
    wellPanel(
      h4("Save leftNotesUp data.frame"),
      helpText(""),
      textInput("leftNotesUpfilename2",
        "File name",
        value = "leftNotesUp",
        width = "50%"
      ),
      uiOutput("buttontableleftNotesUp2")
    )
  }
})


output$searchPanel <- renderUI({
  fluidRow(
    column(
      width = 3,
      br(),
      wellPanel(
        checkboxInput("showWorkflow", "show workflow", value = FALSE)
      ),
      uiOutput("workflowUI"),
      uiOutput("presetUI")
    ),
    column(
      width = 3,
      br(),
      wellPanel(
        checkboxInput("fileInsteadNcbi", "Upload file instead of searching (3a)", value = FALSE)
      ),
      uiOutput("searchUI"),
      uiOutput("chooseFileUI"),
      uiOutput("titleSelectUI"),
      uiOutput("authorsUI")
    ),
    column(
      3,
      br(),
      uiOutput("fetchSelectUI"),
      uiOutput("statsDF")
    ),
    column(
      3,
      br(),
      uiOutput("markColumnUI"),
      uiOutput("loadUI")
    )
  )
})

output$parameterPanel <- renderUI({
  fluidRow(
    column(
      width = 2,
      br(),
      wellPanel(
        style = "padding: 0px 15px 0px 5px",
        fluidRow(
          column(
            4,
            h4("Legend"),
            div(
              title = legendDesc,
              radioButtons("legend", "Pos.", c("aside", "inline", "none"),
                selected = paramValues$legendDefault
              )
            ),
            div(
              style = "margin-top:-12px;margin-bottom:-7px",
              title = remSimiMarkLegDesc,
              checkboxInput(
                "remSimiMarkLeg",
                HTML(paste("remove",
                  "duplicated",
                  sep = "<br/>"
                )),
                paramValues$remSimiMarkLegDefault
              )
            )
          ),
          column(
            4,
            div(
              title = legendWidthDesc,
              numericInput("legendWidth",
                "width",
                paramValues$legendWidthDefault,
                min = 0.25, max = 5, step = 0.05
              )
            ),
            div(
              title = patternDesc,
              textInput(
                "pattern",
                "remove Regex",
                paramValues$patternDefault
              )
            )
          ),
          column(
            4,
            div(
              title = legendHeightDesc,
              numericInput("legendHeight", "height",
                paramValues$legendHeightDefault,
                min = 0.25, max = 5, step = 0.05
              )
            ),
            div(
              title = markNewLineDesc,
              textInput(
                "markNewLine",
                "split sq. marks",
                paramValues$markNewLineDefault
              )
            )
          ) # c
        ) # fR
        , splitLayout(
          cellWidths = c("30%", "2%", "30%", "2%", "30%"),
          div(
            title = forbiddenMarkDesc,
            textInput(
              "forbiddenMark",
              HTML(paste("remove",
                "completely",
                sep = "<br/>"
              )),
              paramValues$forbiddenMarkDefault
            )
          ),
          br(),
          div(
            title = bannedMarkNameDesc,
            textInput(
              "bannedMarkName",
              HTML(paste("ban",
                "this labels",
                sep = "<br/>"
              )),
              paramValues$bannedMarkNameDefault
            )
          ),
          br(),
          div(
            title = bMarkNameAsideDesc,
            checkboxInput(
              "bMarkNameAside",
              tagList(
                div(
                  style = "height:15px;",
                  tags$span("banned", style = "font-size: 80%;")
                ),
                div(
                  style = "height:15px;",
                  tags$span("marks", style = "font-size: 80%;")
                ),
                div(
                  style = "height:15px;",
                  tags$span("to side", style = "font-size: 80%;")
                )
              ),
              paramValues$bMarkNameAsideDefault
            )
          )
        ),
        splitLayout(
          cellWidths = c("30%", "2%", "30%", "2%", "30%"),
          div(
            title = markLabelSizeDesc,
            numericInput("markLabelSize",
              tagList(
                div(
                  style = "height:15px;",
                  tags$span("font")
                ),
                div(
                  style = "height:15px;",
                  tags$span("size")
                )
              ),
              paramValues$markLabelSizeDefault,
              min = 0.1, max = 5, step = 0.05
            )
          ),
          br(),
          div(
            title = markLabelSpacerDesc,
            numericInput("markLabelSpacer",
              tagList(
                div(
                  style = "height:15px;",
                  tags$span("separ.")
                ),
                div(
                  style = "height:15px;",
                  tags$span("from kar.")
                )
              ),
              paramValues$markLabelSpacerDefault,
              min = 0.25, max = 15, step = 0.25
            )
          ),
          br(),
          div(
            title = legendYcoordDesc,
            numericInput("legendYcoord",
              tagList(
                div(
                  style = "height:15px;",
                  tags$span("vertical")
                ),
                div(
                  style = "height:15px;",
                  tags$span("pos.")
                )
              ),
              paramValues$legendYcoordDefault,
              min = 0.0, max = 30, step = 0.2
            )
          )
        ) # sL
      ),
      wellPanel(
        style = "margin-top:-15px;padding: 0px 20px 0px 5px;",
        splitLayout(
          cellWidths = c("28%", "72%"),
          h4("Notes"),
          div(
            title = OTUasNoteDesc,
            checkboxInput(
              "OTUasNote",
              HTML(paste("OTU as note (right)", sep = "<br/>")),
              paramValues$OTUasNoteDefault
            )
          )
        ),
        div(
          style = "margin-top:-10px;",
          title = parseStr2langDesc,
          checkboxInput(
            "parseStr2lang",
            HTML(paste("use str2lang",
              sep = "<br/>"
            )),
            paramValues$parseStr2langDefault
          )
        ),
        splitLayout(
          cellWidths = c("15%", "21%", "21%", "42%"),
          tags$span(),
          tags$span(tags$strong("font")),
          tags$span(),
          tags$span(tags$strong("position"))
        ),
        splitLayout(
          cellWidths = c("15%", "21%", "21%", "21%", "21%"),
          tagList(
            tags$span(tags$strong("right")),
            br(),
            tags$span(tags$strong("notes"))
          ),
          div(
            title = notesTextSizeDesc,
            numericInput("notesTextSize",
              HTML(paste("size", sep = "<br/>")),
              paramValues$notesTextSizeDefault,
              min = 0.1, max = 10, step = 0.1
            )
          ),
          div(
            title = noteFontDesc,
            numericInput("noteFont",
              HTML(paste("style", sep = "<br/>")),
              paramValues$noteFontDefault,
              min = 1, max = 4, step = 1
            )
          ),
          div(
            title = notesPosXDesc,
            numericInput("notesPosX",
              HTML(paste("X-axis", sep = "<br/>")),
              paramValues$notesPosXDefault,
              min = -20, max = 20, step = 0.1
            )
          ),
          div(
            title = notesPosYDesc,
            numericInput("notesPosY",
              HTML(paste("Y-axis", sep = "<br/>")),
              paramValues$notesPosYDefault,
              min = -20, max = 20, step = 0.1
            )
          )
        ),
        splitLayout(
          cellWidths = c("15%", "21%", "21%", "21%", "21%"),
          tagList(
            tags$span(tags$strong("left")),
            br(),
            tags$span(tags$strong("notes"))
          ),
          div(
            title = leftNotesTextSizeDesc,
            numericInput("leftNotesTextSize",
              HTML(paste("size", sep = "<br/>")),
              paramValues$leftNotesTextSizeDefault,
              min = 0.1, max = 10, step = 0.1
            )
          ),
          div(
            title = leftNoteFontDesc,
            numericInput("leftNoteFont",
              HTML(paste("style", sep = "<br/>")),
              paramValues$leftNoteFontDefault,
              min = 1, max = 4, step = 1
            )
          ),
          div(
            title = leftNotesPosXDesc,
            numericInput("leftNotesPosX",
              HTML(paste("X-axis", sep = "<br/>")),
              paramValues$leftNotesPosXDefault,
              min = -5, max = 5, step = 0.1
            )
          ),
          div(
            title = leftNotesPosYDesc,
            numericInput("leftNotesPosY",
              HTML(paste("Y-axis", sep = "<br/>")),
              paramValues$leftNotesPosYDefault,
              min = -5, max = 5, step = 0.1
            )
          )
        ),
        splitLayout(
          cellWidths = c("15%", "21%", "21%", "21%", "21%"),
          tagList(
            tags$span(tags$strong("up left"), style = "font-size: 80%;"),
            br(),
            tags$span(tags$strong("notes"))
          ),
          div(
            title = leftNotesUpTextSizeDesc,
            numericInput("leftNotesUpTextSize",
              HTML(paste("size", sep = "<br/>")),
              paramValues$leftNotesUpTextSizeDefault,
              min = 0.1, max = 10, step = 0.1
            )
          ),
          div(
            title = leftNoteFontUpDesc,
            numericInput("leftNoteFontUp",
              HTML(paste("style", sep = "<br/>")),
              paramValues$leftNoteFontUpDefault,
              min = 1, max = 4, step = 1
            )
          ),
          div(
            title = leftNotesUpPosXDesc,
            numericInput("leftNotesUpPosX",
              HTML(paste("X-axis", sep = "<br/>")),
              paramValues$leftNotesUpPosXDefault,
              min = -5, max = 5, step = 0.1
            )
          ),
          div(
            title = leftNotesUpPosYDesc,
            numericInput("leftNotesUpPosY",
              HTML(paste("Y-axis", sep = "<br/>")),
              paramValues$leftNotesUpPosYDefault,
              min = -5, max = 5, step = 0.1
            )
          )
        )
      ),
      wellPanel(
        style = "margin-top:-15px;padding: 0px 15px 0px 5px",
        h4("Marks"),
        fluidRow(
          column(
            6,
            div(
              style = "margin-top:0px;margin-bottom:-5px",
              title = protrudingDesc,
              numericInput("protruding",
                HTML(paste("protruding of",
                  "cM and arrows",
                  sep = "<br/>"
                )),
                paramValues$protrudingDefault,
                min = 0.05, max = 2.5, step = 0.05
              )
            ),
            div(
              style = "margin-top:0px;margin-bottom:-5px",
              title = arrowheadDesc,
              numericInput("arrowhead",
                HTML(paste("proportion of",
                  "arrowhead",
                  sep = "<br/>"
                )),
                paramValues$arrowheadDefault,
                min = 0.05, max = 1, step = 0.05
              )
            ),
            div(
              title = useOneDotDesc,
              checkboxInput(
                "useOneDot",
                "one dot only",
                paramValues$useOneDotDefault
              )
            ),
            div(
              title = cMBeginCenterDesc,
              checkboxInput(
                "cMBeginCenter",
                HTML(paste("cM mark",
                  "start centered",
                  sep = "<br/>"
                )),
                paramValues$cMBeginCenterDefault
              )
            )
          ),
          column(
            6,
            div(
              style = "margin-top:0px;margin-bottom:-5px",
              title = pMarkFacDesc,
              numericInput("pMarkFac",
                HTML(paste("exProtein",
                  "mark size",
                  sep = "<br/>"
                )),
                paramValues$pMarkFacDefault,
                min = 0.05, max = 1, step = 0.05
              )
            ),
            div(
              title = markDistTypeDesc,
              radioButtons("markDistType", "coord.",
                c(
                  "beginning" = "beg",
                  "center" = "cen"
                ),
                selected = paramValues$markDistTypeDefault
              )
            ),
            div(
              title = originDesc,
              radioButtons("origin", "coord.(hol)",
                c(
                  "bottom" = "b",
                  "top" = "t"
                ),
                selected = paramValues$originDefault
              )
            )
          )
        ),
        div(
          style = "text-align:center;",
          tags$span(tags$strong("thickness"))
        ),
        splitLayout(
          div(
            title = lwd.marksDesc,
            numericInput("lwd.marks",
              HTML(paste("general",
                sep = "<br/>"
              )),
              paramValues$lwd.marksDefault,
              min = 0.25, max = 5, step = 0.25
            )
          ),
          div(
            title = lwd.mimicCenDesc,
            numericInput("lwd.mimicCen",
              HTML(paste("cenStyle",
                sep = "<br/>"
              )),
              paramValues$lwd.mimicCenDefault,
              min = 0.25, max = 5, step = 0.25
            )
          ),
          div(
            title = lwd.cMDesc,
            numericInput("lwd.cM",
              HTML(paste("cM Marks",
                sep = "<br/>"
              )),
              paramValues$lwd.cMDefault,
              min = 0.1, max = 5, step = 0.1
            )
          )
        )
      ),
      wellPanel(
        h4("Add empty kars."),
        fluidRow(
          column(
            8,
            div(
              title = addMissingOTUAfterDesc,
              textInput(
                "addMissingOTUAfter",
                HTML(paste("after OTUs:", sep = "<br/>")),
                paramValues$addMissingOTUAfterDefault
              )
            )
          ),
          column(
            4,
            div(
              title = missOTUspacingsDesc,
              textInput(
                "missOTUspacings",
                HTML(paste("spaces", sep = "<br/>")),
                paramValues$missOTUspacingsDefault
              )
            )
          )
        )
      ),
      wellPanel(
        h4("Special scale"),
        div(
          title = specialOTUNamesDesc,
          textInput(
            "specialOTUNames",
            HTML(paste("OTUs", sep = "<br/>")),
            paramValues$specialOTUNamesDefault
          )
        ),
        splitLayout(
          style = "margin-bottom:-10px;",
          cellWidths = c("33%", "33%", "33%"),
          div(
            title = specialyTitleDesc,
            textInput(
              "specialyTitle",
              HTML(paste("ruler", "title", sep = "<br/>")),
              paramValues$specialyTitleDefault
            )
          ),
          div(
            title = specialChrWidthDesc,
            numericInput("specialChrWidth",
              HTML(paste("chr.", "width", sep = "<br/>")),
              paramValues$specialChrWidthDefault,
              min = 0.1, max = 5, step = 0.1
            )
          ),
          div(
            title = specialChrSpacingDesc,
            numericInput("specialChrSpacing",
              HTML(paste("chr", "spacing", sep = "<br/>")),
              paramValues$specialChrSpacingDefault,
              min = 0.1, max = 5, step = 0.1
            )
          )
        )
      ),
      wellPanel(
        style = "margin-top:-15px;",
        splitLayout(
          cellWidths = c("66%", "33%"),
          div(
            splitLayout(
              cellWidths = c("35%", "60%"),
              h4("OTU"),
              div(
                style = "overflow-x:hidden;",
                title = addOTUNameDesc,
                checkboxInput(
                  "addOTUName",
                  HTML(paste("Show name",
                    sep = "<br/>"
                  )),
                  paramValues$addOTUNameDefault
                )
              )
            ),
            splitLayout(
              cellWidths = c("45%", "45%"),
              div(
                style = "margin-top:20px;overflow-x:hidden;",
                title = OTUfamilyDesc,
                textInput(
                  "OTUfamily",
                  tags$span(tags$strong("Font Family"), style = "font-size: 85%;"),
                  paramValues$OTUfamilyDefault
                )
              ),
              div(
                style = "margin-top:20px;",
                title = OTUTextSizeDesc,
                numericInput("OTUTextSize",
                  HTML(paste("text size",
                    sep = "<br/>"
                  )),
                  paramValues$OTUTextSizeDefault,
                  min = 0.1, max = 5, step = 0.1
                )
              )
            ) # sL
          ),
          div(
            title = OTUfontDesc, style = "font-size: 12px;",
            radioButtons("OTUfont", "font type",
              c("normal" = "1", "bold" = "2", "italics" = "3", "bold italics" = "4"),
              selected = paramValues$OTUfontDefault
            )
          )
        )
      ),
      wellPanel(
        style = "margin-top: -10px",
        splitLayout(
          h4("Layout"),
          div(
            title = circularPlotDesc,
            checkboxInput(
              "circularPlot",
              HTML(paste("Circular", "Plot", sep = "<br/>")),
              paramValues$circularPlotDefault
            )
          )
        ),
        uiOutput("circParam"),
        uiOutput("circParamOTU")
      )
    ),
    column(
      width = 10,
      br(),
      fluidRow(
        uiOutput("colColor"),
        column(
          2,
          wellPanel(
            style = "padding: 0px 5px 0px 5px",
            h4("Chromosomes"),
            fluidRow(
              style = "margin-bottom:-15px;",
              column(
                3,
                div(
                  style = "margin-bottom:-15px",
                  tags$span(tags$strong("order")),
                  tags$span(tags$strong("by:"))
                )
              ),
              column(
                9,
                div(
                  style = "margin-top:-20px;margin-bottom:-15px;",
                  title = orderChrDesc,
                  radioButtons(
                    "orderChr", "",
                    c(
                      "size" = "size",
                      "as in d.f." = "original",
                      "alphab." = "name",
                      "group" = "group",
                      "col. chrNameUp" = "chrNameUp"
                    ),
                    selected = paramValues$orderChrDefault,
                    inline = T
                  )
                )
              )
            ),
            div(
              style = "margin-bottom:-15px;margin-top:-18px;",
              tags$span(tags$strong("Chr. names:"))
            ),
            div(
              style = "margin-top: -10px;",
              title = chrIdDesc,
              radioButtons("chrId",
                "",
                choices = c(
                  "original" = "original", "simple" = "simple",
                  "none" = "none"
                ),
                paramValues$chrIdDefault,
                inline = T
              )
            ),
            splitLayout(
              cellWidths = c("40%", "29%", "29%"),
              div(
                title = chrIdPatternRemDesc,
                textInput(
                  "chrIdPatternRem",
                  HTML(paste(
                    tags$span("remove pattern", style = "font-size: 80%;"),
                    sep = "<br/>"
                  )),
                  paramValues$chrIdPatternRemDefault
                )
              ),
              div(
                title = classChrNameDesc,
                textInput(
                  "classChrName",
                  HTML(paste(
                    tags$span("label"),
                    sep = "<br/>"
                  )),
                  paramValues$classChrNameDefault
                )
              ),
              div(
                title = hideCenLinesDesc,
                numericInput(
                  "hideCenLines",
                  tags$span("hide factor", style = "font-size: 80%;"),
                  paramValues$hideCenLinesDefault
                )
              )
            ),
            splitLayout(
              cellWidths = c("24%", "24%", "24%", "24%"),
              div(
                title = chrWidthDesc,
                numericInput("chrWidth", HTML(paste("", "Width", sep = "</br>")),
                  paramValues$chrWidthDefault,
                  min = 0.1, max = 5, step = 0.05
                )
              ),
              div(
                title = chrSpacingDesc,
                numericInput("chrSpacing", HTML(paste("horiz.", "spacing", sep = "</br>")),
                  paramValues$chrSpacingDefault,
                  min = 0.1, max = 5, step = 0.1
                )
              ),
              div(
                title = groupSeparDesc,
                numericInput("groupSepar", HTML(paste("group", "spacing", sep = "</br>")),
                  paramValues$groupSeparDefault,
                  min = 0.1, max = 5, step = 0.1
                )
              ),
              div(
                title = lwd.chrDesc,
                numericInput("lwd.chr",
                  HTML(paste("border", "width", sep = "</br>")),
                  paramValues$lwd.chrDefault,
                  min = 0, max = 4, step = 0.25
                )
              )
            )
          )
        ),
        column(
          2,
          wellPanel(
            style = "padding: 0px 5px 0px 5px",
            splitLayout(
              cellWidths = c("75%", "24%"),
              div(
                style = "padding: 0px 0px 0px 0px; overflow-x: hidden;",
                h4("Karyotypes"),
                splitLayout(
                  div(
                    title = karSeparDesc,
                    checkboxInput(
                      "karSepar",
                      HTML(paste(
                        tags$span("equally", style = "font-size: 80%;"),
                        tags$span("spaced", style = "font-size: 80%;"),
                        sep = "<br/>"
                      )),
                      paramValues$karSeparDefault
                    )
                  ),
                  div(
                    title = verticalPlotDesc,
                    checkboxInput(
                      "verticalPlot",
                      HTML(paste(
                        tags$span("vertic.", style = "font-size: 80%;"),
                        tags$span("plot", style = "font-size: 80%;"),
                        sep = "<br/>"
                      )),
                      paramValues$verticalPlotDefault
                    )
                  )
                )
              ),
              div(
                title = classGroupNameDesc,
                textInput(
                  "classGroupName",
                  HTML(paste(
                    tags$span("group", style = "font-size: 70%;"),
                    tags$span("label", style = "font-size: 70%;"),
                    sep = "<br/>"
                  )),
                  paramValues$classGroupNameDefault
                )
              )
            ),
            splitLayout(
              style = "margin-bottom:-12px;",
              cellWidths = c("23%", "1%", "23%", "1%", "23%", "1%", "23%"),
              div(
                title = karHeightDesc,
                numericInput("karHeight",
                  tagList(
                    div(
                      style = "height:15px;",
                      tags$span("")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("Height")
                    )
                  ),
                  paramValues$karHeightDefault,
                  min = 0.5, max = 50, step = 0.5
                )
              ),
              br(),
              div(
                title = karHeiSpaceDesc,
                numericInput("karHeiSpace",
                  tagList(
                    div(
                      style = "height:15px;",
                      tags$span("Height", style = "font-size: 70%; ")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("with", style = "font-size: 70%; ")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("space", style = "font-size: 70%; ")
                    )
                  ),
                  paramValues$karHeiSpaceDefault,
                  min = 2, max = 150, step = 0.5
                )
              ),
              br(),
              div(
                title = amoSeparDesc,
                numericInput("amoSepar",
                  tagList(
                    div(
                      style = "height:15px;",
                      tags$span("")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("Vert.")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("separ.")
                    )
                  ),
                  paramValues$amoSeparDefault,
                  min = 0, max = 15, step = 0.5
                )
              ),
              br(),
              div(
                title = karSpaceHorDesc,
                numericInput("karSpaceHor",
                  tagList(
                    div(
                      style = "height:15px;",
                      tags$span("")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("Horiz.")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("separ.")
                    )
                  ),
                  paramValues$karSpaceHorDefault,
                  min = 0, max = 15, step = 0.5
                )
              )
            )
          ),
          wellPanel(
            style = "margin-top:-15px;padding: 0px 5px 0px 5px",
            splitLayout(
              cellWidths = c("50%", "24%", "24%"),
              div(
                h4("Chromatids"),
                splitLayout(
                  cellWidths = c("48%", "48%"),
                  div(
                    title = chromatidsDesc,
                    checkboxInput(
                      "chromatids",
                      tagList(
                        div(
                          style = "height:15px;",
                          tags$span("Show", style = "font-size: 70%;")
                        ),
                        div(
                          style = "height:15px;",
                          tags$span("separ.", style = "font-size: 70%;")
                        )
                      ),
                      paramValues$chromatidsDefault
                    )
                  ),
                  div(
                    title = holocenNotAsChromatidsDesc,
                    checkboxInput(
                      "holocenNotAsChromatids",
                      tagList(
                        div(
                          style = "height:15px;",
                          tags$span("hol.", style = "font-size: 70%;")
                        ),
                        div(
                          style = "height:15px;",
                          tags$span("forbid.", style = "font-size: 65%;")
                        )
                      ),
                      paramValues$holocenNotAsChromatidsDefault
                    )
                  ) # div
                ) # sL
              ) # div
              , div(
                style = "padding: 0px 0px 0px 0px; overflow-x: hidden;",
                title = xModifierDesc,
                numericInput("xModifier",
                  tags$span("Separation", style = "font-size: 80%;"),
                  paramValues$xModifierDefault,
                  min = 2, max = 25, step = 0.5
                )
              ),
              div(
                title = xModMonoHoloRateDesc,
                numericInput("xModMonoHoloRate",
                  tags$span("hol. sep.", style = "font-size: 80%;"),
                  paramValues$xModMonoHoloRateDefault,
                  min = 1, max = 25, step = 1
                )
              )
            ) # sL
          ),
          wellPanel(
            style = "margin-top:-15px; padding: 5px 5px 5px 5px;",
            div(
              class = "innerI",
              splitLayout(
                cellWidths = c("30%", "42%", "30%"),
                actionButton("examDecrease", label = character(0), icon = icon("arrow-circle-left"), class = "btn-success"),
                h4("Example"),
                actionButton("examIncrease", label = character(0), icon = icon("arrow-circle-right"), class = "btn-success")
              )
            )
          )
        ),
        column(
          2,
          wellPanel(
            style = "padding: 0px 5px 0px 5px",
            splitLayout(
              style = "margin-bottom:-10px;",
              cellWidths = c("24%", "24%", "24%", "24%"),
              div(
                style = "overflow-x: hidden;",
                h4("Ruler")
              ),
              div(
                style = "overflow-x: hidden;",
                title = rulerDesc,
                checkboxInput(
                  "ruler",
                  tagList(
                    div(
                      style = "height:15px;",
                      tags$span("Show", style = "font-size: 80%;")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("ruler", style = "font-size: 80%;")
                    )
                  ),
                  paramValues$rulerDefault
                )
              ),
              div(
                style = "overflow-x: hidden;",
                title = useMinorTicksDesc,
                checkboxInput(
                  "useMinorTicks",
                  tagList(
                    div(
                      style = "height:15px;",
                      tags$span("minor", style = "font-size: 80%;")
                    ),
                    div(
                      style = "height:15px;",
                      tags$span("ticks", style = "font-size: 80%;")
                    )
                  ),
                  paramValues$useMinorTicksDefault
                )
              ),
              div(
                style = "overflow-x:hidden;overflow-y:hidden;margin-bottom:0;margin-top:5px;",
                title = miniTickFactorDesc,
                numericInput("miniTickFactor",
                  tagList(
                    div(
                      style = "height:10px;",
                      tags$span("minitick", style = "font-size: 70%;")
                    ),
                    div(
                      style = "height:10px;",
                      tags$span("prop.", style = "font-size: 70%;")
                    )
                  ),
                  paramValues$miniTickFactorDefault,
                  min = 2, max = 20, step = 1
                )
              )
            ),
            splitLayout(
              style = "margin-bottom:0px;margin-top:0px;",
              cellWidths = c("45%", "35%", "15%"),
              span(),
              div( # style="height:15px;",
                span(tags$strong("position", style = "font-size: 90%;")),
                span()
              )
            ),
            splitLayout(
              style = "margin-bottom:-15px;",
              cellWidths = c("22%", "3%", "23%", "23%", "23%"),
              div(
                title = ceilingFactorDesc,
                numericInput("ceilingFactor",
                  tags$span(tags$strong("approx."), style = "font-size: 80%;"),
                  paramValues$ceilingFactorDefault,
                  min = -2, max = 5, step = 1
                )
              ),
              br(),
              div(
                title = rulerPosDesc,
                numericInput("rulerPos", "ruler",
                  paramValues$rulerPosDefault,
                  min = -5, max = 5, step = 0.1
                ),
              ),
              div(
                title = xPosRulerTitleDesc,
                numericInput("xPosRulerTitle",
                  tags$span("title", style = "font-size: 80%;"),
                  paramValues$xPosRulerTitleDefault,
                  min = -5, max = 5, step = 0.1
                )
              ),
              div(
                title = rulerNumberPosDesc,
                numericInput("rulerNumberPos",
                  tags$span("Numbers", style = "font-size: 80%;"),
                  paramValues$rulerNumberPosDefault,
                  min = .1, max = 5, step = 0.1
                )
              )
            ),
            splitLayout(
              style = "margin-bottom:-10px;",
              cellWidths = c("45%", "50%"),
              br(),
              h4(tags$strong("interval"))
            ),
            splitLayout(
              cellWidths = c("22%", "3%", "23%", "23%", "23%"),
              div(
                title = yTitleDesc,
                textInput(
                  "yTitle",
                  HTML(paste("title", sep = "<br/>")),
                  paramValues$yTitleDefault
                )
              ),
              br(),
              div(
                title = rulerIntervalDesc,
                numericInput("rulerInterval",
                  "for \u00B5m",
                  paramValues$rulerIntervalDefault,
                  min = .1, max = 20, step = 0.1
                ),
              ),
              div(
                title = rulerIntervalcMDesc,
                numericInput("rulerIntervalcM",
                  "for cM",
                  paramValues$rulerIntervalcMDefault,
                  min = .1, max = 20, step = 0.1
                ),
              ),
              div(
                title = rulerIntervalMbDesc,
                numericInput("rulerIntervalMb",
                  "for Mb",
                  paramValues$rulerIntervalMbDefault,
                  min = 0.1, max = 10000, step = 0.1
                ),
              )
            ),
            splitLayout(
              style = "margin-bottom:-15px;margin-top:-12px;",
              cellWidths = c("60%", "35%"),
              br(),
              h4(tags$strong("font size"))
            ),
            splitLayout(
              style = "margin-bottom:-10px;",
              cellWidths = c("26%", "24%", "4%", "21%", "21%"),
              div(
                title = ruler.tckDesc,
                numericInput("ruler.tck",
                  tags$span("ticks", style = "font-size: 80%;"),
                  paramValues$ruler.tckDefault,
                  min = -5, max = 5, step = 0.005
                )
              ),
              div(
                title = thresholdDesc,
                numericInput("threshold",
                  tags$span("thresh.", style = "font-size: 80%;"),
                  paramValues$thresholdDefault,
                  min = 5, max = 95, step = 1
                )
              ),
              br(),
              div(
                title = rulerNumberSizeDesc,
                numericInput("rulerNumberSize",
                  tags$span("Numbers", style = "font-size: 70%;"),
                  paramValues$rulerNumberSizeDefault,
                  min = .1, max = 5, step = 0.1
                ),
              ),
              div(
                title = rulerTitleSizeDesc,
                numericInput("rulerTitleSize",
                  "Title",
                  paramValues$rulerTitleSizeDefault,
                  min = .1, max = 5, step = 0.1
                ),
              )
            )
          )
        ),
        column(
          2,
          wellPanel(
            style = "padding: 0px 0px 0px 5px;",
            splitLayout(
              cellWidths = c("45%", "53%"),
              div(
                style = "overflow-x:hidden;",
                h4("Centromere"),
                div(
                  style = "margin-top:-20px",
                  title = cenFormatDesc,
                  radioButtons("cenFormat", "", c("rounded", "triangle", "inProtein"),
                    selected = paramValues$cenFormatDefault
                  )
                )
              ),
              div(
                style = "overflow-x:hidden;overflow-y:hidden;",
                div(
                  title = autoCenSizeDesc,
                  checkboxInput(
                    "autoCenSize",
                    tagList(
                      div( # style="height:15px;",
                        tags$span("auto. size")
                      )
                    ),
                    paramValues$autoCenSizeDefault
                  )
                ),
                div(
                  style = "margin-top:-5px",
                  title = collapseCenDesc,
                  checkboxInput(
                    "collapseCen",
                    tagList(
                      div(
                        tags$span("continuous ruler",
                          style = "font-size: 70%;"
                        )
                      )
                    ),
                    paramValues$collapseCenDefault
                  )
                ),
                splitLayout(
                  style = "margin-bottom:-10px;",
                  cellWidths = c("50%", "48%"),
                  div(
                    title = cenFactorDesc,
                    numericInput("cenFactor",
                      tags$span("modify", style = "font-size: 70%;"),
                      paramValues$cenFactorDefault,
                      min = 0, max = 10, step = 0.25
                    )
                  ),
                  div(
                    title = centromereSizeDesc,
                    numericInput("centromereSize",
                      tags$span("absolute", style = "font-size: 70%;"),
                      paramValues$centromereSizeDefault,
                      min = 0, max = 100000000, step = 1
                    )
                  )
                )
              ) # d
            )
          ),
          wellPanel(
            style = "margin-top:-15px;padding: 0px 20px 0px 5px;",
            splitLayout(
              span(h4("Squareness"), helpText("> 20 = squared")),
              div(
                title = squarenessDesc,
                numericInput("squareness",
                  "",
                  paramValues$squarenessDefault,
                  min = 1, max = 21, step = 0.5
                )
              )
            ),
            splitLayout(
              style = "margin-bottom:-10px;",
              cellWidths = c("50%", "25%", "25%"),
              helpText(
                tags$strong("resolution:"),
                HTML("</br>"),
                tags$span("> = more vertices", style = "font-size: 65%;")
              ),
              div(
                title = markNDesc,
                numericInput("markN",
                  HTML(paste("marks'", sep = "<br/>")),
                  paramValues$markNDefault,
                  min = 1, max = 200, step = 1
                )
              ),
              div(
                title = nDesc,
                numericInput("n",
                  HTML(paste("chr.", sep = "<br/>")),
                  paramValues$nDefault,
                  min = 25, max = 400, step = 1
                )
              )
            )
          ),
          wellPanel(
            style = "margin-top:-15px; padding: 5px 5px 5px 5px;",
            div(
              class = "innerI",
              splitLayout(
                cellWidths = c("30%", "42%", "30%"),
                actionButton("fontDecrease", label = character(0), icon = icon("arrow-alt-circle-down")),
                h4("Fonts' size"),
                actionButton("fontIncrease", label = character(0), icon = icon("arrow-alt-circle-up"))
              )
            )
          )
        ),
        column(
          2,
          wellPanel(
            style = "padding: 0px 5px 0px 5px",
            splitLayout(
              style = "margin-bottom:-10px;",
              cellWidths = c("20%", "38%", "38%"),
              div(
                h4("Plot")
              ),
              div(
                numericInput("widFactor",
                  "width %",
                  80,
                  min = 5, max = 1000, step = 5
                )
              ),
              div(
                numericInput("heiFactor",
                  "height ratio",
                  0.5,
                  min = 0.05, max = 20, step = 0.05
                )
              )
            ),
            splitLayout(
              cellWidths = c("25%", "37%", "38%"),
              br(),
              uiOutput("imageWidth"),
              uiOutput("imageHeight")
            ),
            div(
              style = "text-align:center",
              h5(strong("Format:"))
            ),
            splitLayout(
              style = "margin-top:-7px;",
              radioButtons("pngorsvg", "Display", c("svg", "png"), selected = "svg", inline = T),
              radioButtons("pngorsvgDown", "Download", c("svg", "png"), inline = T)
            ),
            div(
              style = "margin-top:-7px;margin-bottom:7px",
              downloadButton("pngorsvgDownButton",
                tags$strong("Download plot"),
                class = "btn-info"
              )
            )
          ) # wP
          , wellPanel(
            style = "margin-top:-15px; padding: 2px 5px 0px 5px",
            helpText(presetUseText),
            div(
              style = "margin-top:-4px;margin-bottom:5px;",
              uiOutput("buttonPresets2")
            )
          )
        ) # c
      ),
      fluidRow(
        column(
          10,
          wellPanel(
            style = "margin-top:-15px;overflow-y:auto;text-align:center;",
            div(
              style = "overflow-y:auto;overflow-x:auto;min-width:1030px;",
              imageOutput("idiogramPlot",
                height = "auto"
              )
            )
          )
        ),
        uiOutput("lastCol")
      )
    )
  )
})


output$lastCol <- renderUI({
  column(
    width = 2,
    wellPanel(
      style = "margin-top:-15px;padding: 5px 5px 5px 5px",
      div(
        style = "margin-bottom:-10px",
        splitLayout(
          cellWidths = c("30%", "20%", "30%", "20%"),
          h5(tags$strong("Zoom (X)")),
          actionButton("left", label = character(0), icon = icon("search-minus")),
          numericInput("hwModifier", label = NULL, min = 0, max = 10, value = 1, step = 0.1),
          actionButton("right", label = character(0), icon = icon("search-plus"))
        )
      )
    ),
    wellPanel(
      style = "margin-top:-15px;padding: 0px 15px 0px 5px",
      h4("Margins"),
      splitLayout(
        style = "margin-bottom:-10px;",
        cellWidths = c("25%", "25%", "25%", "25%"),
        div(
          title = ylimBotModDesc,
          numericInput("ylimBotMod", "bottom",
            paramValues$ylimBotModDefault,
            min = -5, max = 5, step = 0.5
          )
        ),
        div(
          title = ylimTopModDesc,
          numericInput("ylimTopMod", "top",
            paramValues$ylimTopModDefault,
            min = -5, max = 5, step = 0.5
          )
        ),
        div(
          title = xlimLeftModDesc,
          numericInput("xlimLeftMod", "left",
            paramValues$xlimLeftModDefault,
            min = -5, max = 5, step = 0.5
          )
        ),
        div(
          title = xlimRightModDesc,
          numericInput("xlimRightMod", "right",
            paramValues$xlimRightModDefault,
            min = -6, max = 6, step = 0.4
          )
        )
      )
    ),
    wellPanel(
      style = "margin-top:-15px;padding: 0px 15px 0px 5px",
      h5("Move all karyotypes"),
      splitLayout(
        div(
          title = moveAllKarValueHorDesc,
          numericInput("moveAllKarValueHor",
            HTML(paste("horiz.", sep = "<br/>")),
            paramValues$moveAllKarValueHorDefault,
            min = 0.05, max = 5, step = 0.05
          )
        ),
        div(
          title = moveAllKarValueYDesc,
          numericInput("moveAllKarValueY",
            HTML(paste("vertical", sep = "<br/>")),
            paramValues$moveAllKarValueYDefault,
            min = 0.05, max = 5, step = 0.05
          )
        )
      )
    ),
    wellPanel(
      style = "margin-top:-15px;padding: 0px 15px 0px 5px",
      h4("Indices & Info."),
      splitLayout(
        div(
          title = indexIdTextSizeDesc,
          numericInput("indexIdTextSize",
            HTML(paste("font", "size", sep = "<br/>")),
            paramValues$indexIdTextSizeDefault,
            min = 0.05, max = 5, step = 0.05
          )
        ),
        div(
          title = distTextChrDesc,
          numericInput("distTextChr",
            HTML(paste(
              tags$span(tags$strong("chr. to text"), style = "font-size: 70%;"),
              tags$span(tags$strong("separation"), style = "font-size: 70%;"),
              sep = "<br/>"
            )),
            paramValues$distTextChrDefault,
            min = 0.2, max = 5, step = 0.1
          )
        ),
        div(
          style = "overflow-x:hidden",
          title = karIndexDesc,
          checkboxInput(
            "karIndex",
            tagList(
              div(
                style = "height:15px;",
                tags$span("Add kar.", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("ind. A/A2", style = "font-size: 70%;")
              )
            ),
            paramValues$karIndexDefault
          )
        ) # d
      ),
      splitLayout(
        div(
          title = morphoDesc, style = "font-size: 13px;",
          radioButtons("morpho", "Morphology",
            c(
              "Guerra (1986)" = "Guerra",
              "Levan (1974)" = "Levan",
              "both" = "both",
              "none" = "none"
            ),
            selected = paramValues$morphoDefault
          )
        ),
        div(
          title = chrIndexDesc, style = "font-size: 13px;",
          radioButtons("chrIndex", "Add AR & CI",
            c(
              "Chr. Index" = "CI",
              "Arm Ratio" = "AR",
              "both" = "both",
              "none" = "none"
            ),
            selected = paramValues$chrIndexDefault
          )
        )
      ),
      splitLayout(
        style = "margin-bottom:-10px;",
        cellWidths = c("33%", "67%"),
        br(),
        h4(tags$strong("Show"))
      ),
      splitLayout(
        div(
          title = chrSizeDesc,
          checkboxInput(
            "chrSize",
            tagList(
              div(
                style = "height:15px;",
                tags$span("Chrom.", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("size", style = "font-size: 70%;")
              )
            ),
            paramValues$chrSizeDefault
          )
        ),
        div(
          style = "overflow-x:hidden;",
          title = chrNameUpDesc,
          checkboxInput(
            "chrNameUp",
            tagList(
              div(
                style = "height:15px;",
                tags$span("Chrom.", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("name up", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$p("(column", style = "font-size: 80%;")
              ),
              div(
                style = "height:15px;",
                tags$p("required)", style = "font-size: 80%;")
              )
            ),
            paramValues$chrNameUpDefault
          )
        ),
        div(
          title = chrSizeMbpDesc,
          checkboxInput(
            "chrSizeMbp",
            tagList(
              div(
                style = "height:15px;",
                tags$span("Chr.")
              ),
              div(
                style = "height:15px;",
                tags$span("size")
              ),
              div(
                style = "height:15px;",
                tags$span("in Mbp")
              ),
              div(
                style = "height:15px;",
                tags$p("(column", style = "font-size: 80%;")
              ),
              div(
                style = "height:15px;",
                tags$p("required)", style = "font-size: 80%;")
              )
            ),
            paramValues$chrSizeMbpDefault
          )
        )
      ) # sL
      , splitLayout(
        cellWidths = c("24%", "24%", "24%", "24%"),
        div(
          title = classChrNameUpDesc,
          textInput(
            "classChrNameUp",
            tagList(
              div(
                style = "height:15px;",
                tags$span("up", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("label", style = "font-size: 70%;")
              )
            ),
            paramValues$classChrNameUpDefault
          )
        ),
        div(
          style = "overflow-x:hidden",
          title = nsmallDesc,
          numericInput(
            "nsmall",
            tagList(
              div(
                style = "height:15px;",
                tags$span("approx.", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("values", style = "font-size: 70%;")
              )
            ),
            paramValues$nsmallDefault
          )
        ),
        div(
          title = karIndexPosDesc,
          numericInput(
            "karIndexPos",
            tagList(
              div(
                style = "height:15px;",
                tags$span("kar.", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("ind. pos.", style = "font-size: 70%;")
              )
            ),
            paramValues$karIndexPosDefault
          )
        ),
        div(
          title = nameChrIndexPosDesc,
          numericInput(
            "nameChrIndexPos",
            tagList(
              div(
                style = "height:15px;",
                tags$span("chr.", style = "font-size: 70%;")
              ),
              div(
                style = "height:15px;",
                tags$span("ind. pos.", style = "font-size: 70%;")
              )
            ),
            paramValues$nameChrIndexPosDefault
          )
        )
      ),
      splitLayout(
        div(
          title = perAsFractionDesc,
          checkboxInput(
            "perAsFraction",
            HTML(paste("% as", "fraction", sep = "<br/>")),
            paramValues$perAsFractionDefault
          )
        ),
        div(
          title = showMarkPosDesc,
          checkboxInput(
            "showMarkPos",
            HTML(paste("mark's", "pos.", sep = "<br/>")),
            paramValues$showMarkPosDefault
          )
        )
      ),
      splitLayout(
        div(
          title = markPerDesc,
          textInput(
            "markPer",
            HTML(paste("% of this", "marks:", sep = "<br/>")),
            paramValues$markPerDefault
          )
        ),
        div(
          title = bToRemoveDesc,
          textInput(
            "bToRemove",
            HTML(paste("hide pos. of", "this marks", sep = "<br/>")),
            paramValues$bToRemoveDefault
          )
        )
      )
    ),
    wellPanel(
      style = "margin-top:-15px;padding: 0px 20px 0px 5px",
      splitLayout(
        cellWidths = c("30%", "34%"),
        h4("Anchor"),
        div(
          title = anchorDesc,
          checkboxInput(
            "anchor",
            HTML(paste("Show", "anchor", sep = "<br/>")),
            paramValues$anchorDefault
          )
        )
      ),
      splitLayout(
        div(
          title = anchorVsizeFDesc,
          numericInput("anchorVsizeF",
            HTML(paste("Mod.", "anchor-y", sep = "<br/>")),
            paramValues$anchorVsizeFDefault,
            min = .2, max = 3, step = 0.2
          )
        ),
        div(
          title = anchorHsizeFDesc,
          numericInput("anchorHsizeF",
            HTML(paste("Mod.", "anchor-x", sep = "<br/>")),
            paramValues$anchorHsizeFDefault,
            min = .2, max = 3, step = 0.2
          )
        ),
        div(
          title = mkhValueDesc,
          numericInput("mkhValue",
            HTML(paste("move", "kar.", sep = "<br/>")),
            paramValues$mkhValueDefault,
            min = -10, max = 10, step = 0.25
          )
        )
      ),
      splitLayout(
        div(
          title = karAnchorLeftDesc,
          textInput("karAnchorLeft",
            HTML(paste("kar.", "to anchor", sep = "<br/>")),
            value  = paramValues$karAnchorLeftDefault,
            width = "100%"
          )
        ),
        div(
          title = moveKarHorDesc,
          textInput("moveKarHor",
            HTML(paste("kar.", "to move", sep = "<br/>")),
            value  = paramValues$moveKarHorDefault,
            width = "100%"
          )
        )
      ),
      splitLayout(
        div(
          title = moveAnchorVDesc,
          numericInput("moveAnchorV", "move ver.", paramValues$moveAnchorVDefault, min = -10, max = 10, step = 0.25)
        ),
        div(
          title = moveAnchorHDesc,
          numericInput("moveAnchorH", "move hor.", paramValues$moveAnchorHDefault, min = -10, max = 10, step = 0.25)
        )
      )
    )
  )
})

output$chrColorUI <- renderUI({
  tagList(
    splitLayout(
      cellWidths = c("50%", "50%"),
      h4("Colors"),
      h5("border")
    ),
    splitLayout(
      style = "margin-top:-10px;margin-bottom:-10px",
      cellWidths = c("24%", "24%", "24%", "24%"),
      div(
        title = cenColorDesc,
        textInput("cenColor", "centr.", paramValues$cenColorDefault)
      ),
      div(
        title = chrColorDesc,
        textInput("chrColor", "chr.", paramValues$chrColorDefault)
      ),
      div(
        title = colorBorderMarkDesc,
        textInput(
          "colorBorderMark",
          HTML("Marks"),
          paramValues$colorBorderMarkDefault
        )
      ),
      div(
        title = chrBorderColorDesc,
        textInput("chrBorderColor", "chr.", paramValues$chrBorderColorDefault)
      )
    ),
    splitLayout(
      div(
        title = fixCenBorderDesc,
        checkboxInput(
          "fixCenBorder",
          tagList(
            tags$p("chr. color to", style = "font-size: 70%; font-weight: normal;"),
            tags$p("cen. border color", style = "font-size: 70%; font-weight: normal;")
          ),
          paramValues$fixCenBorderDefault
        )
      ),
      div(
        title = gishCenBorderDesc,
        checkboxInput(
          "gishCenBorder",
          HTML(paste("mark color", "to border", sep = "<br/>")),
          paramValues$gishCenBorderDefault
        )
      )
    )
  )
})

output$marksUIcolor <- renderUI({
  tagList(
    splitLayout(
      cellWidths = c("70%", "29%"),
      div(
        style = "padding:0px 2px 0px 2px; margin-top:0px;",
        title = mycolorsDesc,
        textInput(
          "mycolors",
          HTML(paste(
            "Marks",
            tags$p("(optional, see data.frame page)", style = "font-size: 70%; font-weight: normal;"),
            tags$p("here, only comma separated", style = "font-size: 80%;")
          )),
          paramValues$mycolorsDefault
        )
      ),
      div(
        style = "padding:0px 2px 2px 2px; margin-top: 41px;",
        title = alpha_valDesc,
        numericInput("alpha_val",
          HTML(paste(tags$p("transparency", style = "font-size: 70%;"))),
          paramValues$alpha_valDefault,
          min = 0, max = 1, step = 0.1
        )
      )
    )
  )
})

output$colColor <- renderUI({
  column(
    width = 2,
    wellPanel(
      style = "background: #F7DCDA; padding: 0px 5px 0px 5px",
      uiOutput("chrColorUI"),
      uiOutput("marksUIcolor")
    ),
    wellPanel(
      style = "margin-top:-15px;padding: 8px 5px 9px 5px",
      textOutput(outputId = "currentExample")
    )
  )
})
