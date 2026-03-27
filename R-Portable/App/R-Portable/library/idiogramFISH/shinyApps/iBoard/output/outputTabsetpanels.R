
output$tabsetpanel0UI <- renderUI({
  tabsetPanel(
    id = "tabsetpanel0",
    tabPanel("1. Load example",
      tags$head(tags$style(HTML(" wpanel .tab-content {margin-left:50px;}"))),
      value = "exampleTab",
      div(
        class = "wpanel",
        uiOutput("examplepanel")
      )
    )
  )
})

output$tabsetpanel1UI <- renderUI({
  tabsetPanel(
    id = "tabsetpanel1",
    tabPanel("1. Chr. data data.frame",
      value = "dfChrTab",
      uiOutput("dfchrpanel")
    ),
    tabPanel("2. Marks' pos. d.f.",
      value = "dfMarkTab",
      uiOutput("dfmarkpanel")
    ),
    tabPanel("3. Marks' style d.f.",
      value = "dfMSTab",
      uiOutput("dfMStylepanel")
    ),
    tabPanel("4. Notes' data.frames",
      value = "notesTab",
      uiOutput("dfnotespanel")
    )
  )
})

output$tabsetpanel5UI <- renderUI({
  tabsetPanel(
    id = "tabsetpanel5",
    tabPanel("Indices",
      value = "indicesTab",
      uiOutput("indicespanel")
    ),
    tabPanel("Marks",
      value = "marksTab",
      uiOutput("markspanel")
    )
  )
})

output$tabsetpanel2UI <- renderUI({
  tabsetPanel(
    id = "tabsetpanel2",
    tabPanel("1. Parameters & Plot",
      value = "paramTab"
      # ,style= "min-width:1366px;max-width:1920px;overflow:auto"
      , uiOutput("parameterPanel")
    ),
    tabPanel("2. Log",
      value = "logTab",
      uiOutput("logpanel")
    ),
    tabPanel("3. code",
      value = "codeTab",
      uiOutput("strpanel")
    )
  )
})

output$tabsetpanel4UI <- renderUI({
  tabsetPanel(
    id = "tabsetpanel4",
    tabPanel("1. Search",
      value = "searchTab",
      uiOutput("searchPanel")
    )
  )
})
