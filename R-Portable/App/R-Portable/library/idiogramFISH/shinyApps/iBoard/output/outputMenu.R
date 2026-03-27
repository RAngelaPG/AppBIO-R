
#
#   left submenus
#

output$mysidebar <- renderUI({
  sidebarMenu(
    id = "tabs",
    sidebarMenuOutput("menu1")
  )
})

output$menu1 <- renderMenu({
  sidebarMenu(
    menuItem("Examples (presets)", tabName = "examplesMenu", icon = icon("align-left")),
    menuItem("Nucleotides", tabName = "nuccoreMenu", icon = icon("dna")),
    menuItem("data.frames (input)", tabName = "DFsMenu", icon = icon("table")),
    menuItem("Parameters & Plot",
      tabName = "parameterPlotMenu", icon = icon("chalkboard"),
      selected = TRUE
    ),
    menuItem("stats (output)", tabName = "indicesMenu", icon = icon("table")),
    menuItem("About", tabName = "aboutMenu", icon = icon("scroll"))
  )
}) # output
