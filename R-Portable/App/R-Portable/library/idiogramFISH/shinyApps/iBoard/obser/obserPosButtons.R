#
#   update currentM and current tab reactive values
#

observeEvent(
  input[["tabsetpanel1"]],
  {
    Current$Tab <- input[["tabsetpanel1"]]
  }
)

observeEvent(
  input[["tabsetpanel2"]],
  {
    Current$Tab2 <- input[["tabsetpanel2"]]
  }
)

observeEvent(
  input[["tabsetpanel4"]],
  {
    Current$Tab4 <- input[["tabsetpanel4"]]
  }
)

observeEvent(
  input[["tabsetpanel5"]],
  {
    Current$Tab5 <- input[["tabsetpanel5"]]
  }
)

observeEvent(
  input[["tabs"]],
  {
    CurrentM$menu <- input[["tabs"]]
  }
)


#
#   jump among submenus or tabs
#

#
observeEvent(
  {
    input[["jumpToPrevMenu0"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) - 1
    if (tab_id_position == 0) tab_id_position <- length(menulist)
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

observeEvent(
  {
    input[["jumpToPrevMenu"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) - 1
    if (tab_id_position == 0) tab_id_position <- length(menulist)
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)
observeEvent(
  {
    input[["jumpToPrevMenu2"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) - 1
    if (tab_id_position == 0) tab_id_position <- length(menulist)
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)
observeEvent(
  {
    input[["jumpToPrevMenu3"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) - 1
    if (tab_id_position == 0) tab_id_position <- length(menulist)
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

observeEvent(
  {
    input[["jumpToPrevMenu4"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) - 1
    if (tab_id_position == 0) tab_id_position <- length(menulist)
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

observeEvent(
  {
    input[["jumpToPrevMenu5"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) - 1
    if (tab_id_position == 0) tab_id_position <- length(menulist)
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

#
#   next
#

observeEvent(
  input[["jumpToNextMenu0"]],
  {
    tab_id_position <- match(CurrentM$menu, menulist) + 1
    if (tab_id_position > length(menulist)) tab_id_position <- 1
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1]) # avoids log error
  }
)

observeEvent(
  input[["jumpToNextMenu"]],
  {
    tab_id_position <- match(CurrentM$menu, menulist) + 1
    if (tab_id_position > length(menulist)) tab_id_position <- 1
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

observeEvent(
  input[["jumpToNextMenu2"]],
  {
    tab_id_position <- match(CurrentM$menu, menulist) + 1
    if (tab_id_position > length(menulist)) tab_id_position <- 1
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

observeEvent(
  {
    input[["jumpToNextMenu3"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) + 1
    if (tab_id_position > length(menulist)) tab_id_position <- 1
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1]) # avoids log error
  }
)

observeEvent(
  {
    input[["jumpToNextMenu4"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) + 1
    if (tab_id_position > length(menulist)) tab_id_position <- 1
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1])
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

observeEvent(
  {
    input[["jumpToNextMenu5"]]
  },
  {
    tab_id_position <- match(CurrentM$menu, menulist) + 1
    if (tab_id_position > length(menulist)) tab_id_position <- 1
    CurrentM$menu <- menulist[tab_id_position]
    updateTabItems(session, "tabs", menulist[tab_id_position])
    updateTabItems(session, "tabsetpanel2", tablist2[1])
    updateTabItems(session, "tabsetpanel5", tablist5[1])
  }
)

#
#   tab jumps
#

observeEvent(
  input[["jumpToNext"]],
  {
    tab_id_position <- match(Current$Tab, tablist) + 1
    if (tab_id_position > length(tablist)) tab_id_position <- 1
    Current$Tab <- tablist[tab_id_position]
    updateTabItems(session, "tabsetpanel1", tablist[tab_id_position])
  }
)

observeEvent(
  input[["jumpToPrev"]],
  {
    tab_id_position <- match(Current$Tab, tablist) - 1
    if (tab_id_position == 0) tab_id_position <- length(tablist)
    Current$Tab <- tablist[tab_id_position]
    updateTabItems(session, "tabsetpanel1", tablist[tab_id_position])
  }
)

observeEvent(
  input[["jumpToNext2"]],
  {
    tab_id_position <- match(Current$Tab2, tablist2) + 1
    if (tab_id_position > length(tablist2)) tab_id_position <- 1
    Current$Tab2 <- tablist2[tab_id_position]
    updateTabItems(session, "tabsetpanel2", tablist2[tab_id_position])
  }
)

observeEvent(
  input[["jumpToPrev2"]],
  {
    tab_id_position <- match(Current$Tab2, tablist2) - 1
    if (tab_id_position == 0) tab_id_position <- length(tablist2)
    Current$Tab2 <- tablist2[tab_id_position]
    updateTabItems(session, "tabsetpanel2", tablist2[tab_id_position])
  }
)

observeEvent(
  input[["jumpToNext5"]],
  {
    tab_id_position <- match(Current$Tab5, tablist5) + 1
    if (tab_id_position > length(tablist5)) tab_id_position <- 1
    Current$Tab5 <- tablist5[tab_id_position]
    updateTabItems(session, "tabsetpanel5", tablist5[tab_id_position])
  }
)

observeEvent(
  input[["jumpToPrev5"]],
  {
    tab_id_position <- match(Current$Tab5, tablist5) - 1
    if (tab_id_position == 0) tab_id_position <- length(tablist5)
    Current$Tab5 <- tablist5[tab_id_position]
    updateTabItems(session, "tabsetpanel5", tablist5[tab_id_position])
  }
)

observeEvent(input$loadDFbutton, {
  Sys.sleep(0.1)
  CurrentM$menu <- menulist[4] # plot
  updateTabItems(session, "tabs", menulist[4])
  updateTabItems(session, "tabsetpanel2", tablist2[1]) # paramTab
})

observeEvent(input$exampleButton, {
  Sys.sleep(0.3)

  sel <- input$exampleId
  sel <- gsub(".*-", "", sel)

  if (as.numeric(sel) %in% 9:10) {
    CurrentM$menu <- menulist[2]
    updateTabItems(session, "tabs", menulist[2])
    updateTabItems(session, "tabsetpanel4", tablist4[1])
  } else {
    CurrentM$menu <- menulist[4]
    updateTabItems(session, "tabs", menulist[4])
    updateTabItems(session, "tabsetpanel2", tablist2[1])
  }
})

observeEvent(input$upPresetButton, {
  Sys.sleep(0.3)
  CurrentM$menu <- menulist[4]
  updateTabItems(session, "tabs", menulist[4])
  updateTabItems(session, "tabsetpanel2", tablist2[1])
})



# input$examIncrease
# input$examDecrease
