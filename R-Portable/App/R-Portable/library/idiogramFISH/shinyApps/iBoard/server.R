options(shiny.maxRequestSize = 30 * 1024^2)

tablist0 <- c("exampleTab")

tablist <- c("dfChrTab", "dfMarkTab", "dfMSTab", "notesTab")

tablist2 <- c("paramTab", "logTab", "codeTab")

tablist4 <- c("searchTab")

tablist5 <- c("indicesTab", "marksTab")

menulist <- c("examplesMenu",
  "nuccoreMenu",
  "DFsMenu",
  "parameterPlotMenu",
  "indicesMenu",
  "aboutMenu"
)

add_row <- function(chapter_note, chapter_name) {
  tags$tr(
    tags$td(
      tags$a(href = paste0("https://colab.research.google.com/github/fernandoroa/IFjupyter/blob/main/", chapter_note),
        HTML(paste0(chapter_name, "&emsp;&emsp;")),
        target = "_blank")
    ),
    tags$td(
      tags$a(href = paste0("https://github.com/fernandoroa/IFjupyter/blob/main/", chapter_note),
        tags$kbd("link"),
        target = "_blank"
      )
    ),
    tags$td(
      tags$a(href = paste0("https://github.com/fernandoroa/IFjupyter/raw/main/", chapter_note),
        tags$kbd("Raw"),
        target = "_blank"
      )
    )
  )
}

server <- function(input, output, session) {

  print(paste0("Running in: ", isolate(session$clientData$url_hostname), ":", isolate(session$clientData$url_port)))

  #
  # values

  source(file.path("presetVals.R"),         local = TRUE, chdir = TRUE)

  source(file.path("reactiveV.R"),          local = TRUE, chdir = TRUE)

  #
  # outputs

  source(file.path("output/outputMenu.R"),         local = TRUE, chdir = TRUE)

  source(file.path("output/outputTabsetpanels.R"), local = TRUE, chdir = TRUE)

  source(file.path("output/outputTabs.R"),         local = TRUE, chdir = TRUE)

  source(file.path("output/outputPanels.R"),       local = TRUE, chdir = TRUE)

  source(file.path("output/outputNuccore.R"),      local = TRUE, chdir = TRUE)

  source(file.path("output/outputDataframes.R"),   local = TRUE, chdir = TRUE)

  source(file.path("output/outputPlot.R"),         local = TRUE, chdir = TRUE)

  source(file.path("output/outputLogCode.R"),      local = TRUE, chdir = TRUE)

  #
  # observers
  source(file.path("obser/obserPosButtons.R"),    local = TRUE, chdir = TRUE)

  source(file.path("obser/obserPresets.R"),       local = TRUE, chdir = TRUE)

  source(file.path("obser/obserNuccore.R"),       local = TRUE, chdir = TRUE)

  source(file.path("obser/obserDataframes.R"),    local = TRUE, chdir = TRUE)

  source(file.path("obser/obserSwap.R"),          local = TRUE, chdir = TRUE)

  source(file.path("obser/obserColumnMod.R"),     local = TRUE, chdir = TRUE)

  source(file.path("obser/obserPlot.R"),          local = TRUE, chdir = TRUE)

  source(file.path("obser/obserIndices.R"),          local = TRUE, chdir = TRUE)


  outputOptions(output, "tabsetpanel1UI", suspendWhenHidden = FALSE)
  outputOptions(output, "tabsetpanel2UI", suspendWhenHidden = FALSE)
  outputOptions(output, "tabsetpanel4UI", suspendWhenHidden = FALSE)

  outputOptions(output, "myTabs", suspendWhenHidden = FALSE)

  outputOptions(output, "presetUI",       suspendWhenHidden = FALSE)
  outputOptions(output, "examplepanel",   suspendWhenHidden = FALSE)

  outputOptions(output, "dfchrpanel",     suspendWhenHidden = FALSE)
  outputOptions(output, "out",            suspendWhenHidden = FALSE)

  outputOptions(output, "dfmarkpanel",    suspendWhenHidden = FALSE)
  outputOptions(output, "outMark",        suspendWhenHidden = FALSE)

  outputOptions(output, "dfMStylepanel",  suspendWhenHidden = FALSE)
  outputOptions(output, "outMStyle",      suspendWhenHidden = FALSE)

  outputOptions(output, "dfnotespanel",   suspendWhenHidden = FALSE)
  outputOptions(output, "outNotes",       suspendWhenHidden = FALSE)

  outputOptions(output, "parameterPanel", suspendWhenHidden = FALSE)
  outputOptions(output, "logpanel",       suspendWhenHidden = FALSE)
  outputOptions(output, "strpanel",       suspendWhenHidden = FALSE)

  outputOptions(output, "searchPanel",    suspendWhenHidden = FALSE)
  outputOptions(output, "searchUI",       suspendWhenHidden = FALSE)
  outputOptions(output, "markColumnUI",   suspendWhenHidden = FALSE)
  outputOptions(output, "authorsUI",      suspendWhenHidden = FALSE)

  outputOptions(output, "saveleftNotesUI",      suspendWhenHidden = FALSE)
  outputOptions(output, "saveleftNotesUpUI",      suspendWhenHidden = FALSE)
  outputOptions(output, "saveNotesUI",      suspendWhenHidden = FALSE)

  outputOptions(output, "saveChrSizeUI",      suspendWhenHidden = FALSE)
  outputOptions(output, "saveMarkPosUI",      suspendWhenHidden = FALSE)
  outputOptions(output, "saveMarkStyleUI",      suspendWhenHidden = FALSE)


}
