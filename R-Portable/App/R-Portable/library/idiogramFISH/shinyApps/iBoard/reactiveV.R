
values <- reactiveValues(number = 0,
  mylistAllClass = "numeric",
  ARCImessage = "",
  df1 = data.frame(),
  df1Mark = data.frame(),
  ARCI = data.frame(),
  notes = data.frame(),
  leftNotes = data.frame(),
  leftNotesUp = data.frame(),
  termButtonVal = 0,
  entrez_titles = NA,
  entrez_selected = NA,
  counter = 0,
  decision = "",
  state = "ready",
  pngorsvg = "svg",
  plot = TRUE,
  stop = FALSE,
  rentrezPkg = TRUE,
  renInstall = FALSE,
  errorMessage = "Search failed, change string or check internet",
  renMiss = "unable, rentrez package missing",
  current_len = 0,
  paramVec = paramVec,
  canButton = TRUE,
  maxEx = maxEx,
  go = TRUE
)

#
#   first tab

Current <- reactiveValues(Tab = "dfChrTab", Tab2 = "paramTab", Tab4 = "searchTab", Tab5 = "indicesTab")

#
# first submenu, see local

CurrentM <- reactiveValues(
  menu = "parameterPlotMenu"
)

scriptR <- reactive({
  df <- values[["strFun"]]
})

presetsR <- reactive({
  list1 <- values[["presets"]]
})

current_lenR <- reactive({
  v <- values[["current_len"]]
})

if (suppressWarnings(system.file(package = "rentrez") == "")) {
  values[["rentrezPkg"]] <- FALSE
}

numberRe <- reactive({
  n <- values[["number"]]
})
