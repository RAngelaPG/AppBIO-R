# nucFile updates reactives
# termButton  updates reactives


observeEvent(input$nucFile, {
  showModal(modalDialog(
    title = "3a. Reading file, please wait",
    "this pop-up will close after completion. Press ESC to wait in shiny app",
    easyClose = TRUE,
    footer = modalButton("Wait in shiny app")
  ))

  mylistAll <- 0
  values[["mylistAllClass"]] <- class(mylistAll)

  values[["termButtonVal"]] <- 0
  values[["searchStatus"]] <- FALSE
  values[["button3ab"]] <- NA

  inFile <- input[["nucFile"]]

  validate(
    need(try({
      mylistAll <- genBankReadIF(inFile$datapath)
      class(mylistAll) == "list"
    }), "reading")
  )

  values[["mylistAllClass"]] <- class(mylistAll)

  remove(inFile)

  idSearch <- mylistAll$gbdfMain$value[which(mylistAll$gbdfMain$field == "ACCESSION")]

  values[["authors"]] <- p(paste0(idSearch, ": ", mylistAll$gbdfMain[which(mylistAll$gbdfMain$field == "AUTHORS"), ][1, 2]))

  #
  #   generate d.f. chr. data
  #

  chrDF <- data.frame(chrName = 1, chrSize = mylistAll$source$end)

  mylistAll$source <- NULL

  values[["fetch_listAll"]] <- mylistAll

  chrDF$OTU <- mylistAll$gbdfMain[which(mylistAll$gbdfMain$field == "DEFINITION"), ]$value

  values[["geneChrDF"]] <- chrDF

  genic <- grep("gbdf|source", names(mylistAll), value = TRUE, invert = TRUE)

  if (length(genic) == 0) {
    values[["annotated"]] <- FALSE
    removeModal()
    return()
  } else {
    values[["annotated"]] <- TRUE
  }

  mylist <- mylistAll[which(names(mylistAll) %in% genic)]

  values[["names_fetch_list"]] <- names(mylist)
  values[["fetch_list_selected"]] <- grep("CDS", names(mylist), invert = TRUE)
  values[["fetch_list"]] <- mylist

  validate(
    need(try(!invalid(values[["names_fetch_list"]])), "downloading 912")
  )
  values[["button3ab"]] <- 10
  values[["termButtonVal"]] <- 1
  values[["searchStatus"]] <- TRUE
  gc()
  removeModal()
})

##############################################################################
#
#   search nuccore
#
##############################################################################

observeEvent(input$termButton, {
  if (values[["rentrezPkg"]] == FALSE) {
    showModal(modalDialog(
      title = "For this to work, rentrez needs to be installed",
      tagList("Do you want to install rentrez?"),
      easyClose = FALSE,
      footer = list(
        actionButton("installRen", "Yes, install"),
        actionButton("dontInst", "No, I will not search")
      )
    ))
  } else {
    values[["renInstall"]] <- FALSE
    values[["errorMessage"]] <- "Search failed, change string or check internet"
  }
})

observeEvent(input$installRen, {
  if (system.file(package = "rentrez") == "") {
    tryCatch(detach("package:rentrez", unload = TRUE), error = function(e) {
      ""
    })
    if (!"rentrez" %in% (.packages())) {
      message("installing rentrez")
      tryCatch(utils::install.packages("rentrez"), error = function(w) {
        message("failure installing rentrez")
        return("failure")
      })
    } else {
      message("rentrez loaded, aborting")
    }
  }

  if (system.file(package = "rentrez") == "") {
    values[["rentrezPkg"]] <- FALSE
    values[["renInstall"]] <- FALSE
    values[["renMiss"]] <- "check Internet"
    values[["errorMessage"]] <- "check internet"
  } else {
    values[["rentrezPkg"]] <- TRUE
    values[["renInstall"]] <- TRUE
    values[["errorMessage"]] <- "Try again, press 2."
  }
  removeModal()
})

observeEvent(input$dontInst, {
  if (system.file(package = "rentrez") == "") {
    values[["rentrezPkg"]] <- FALSE
    values[["renMiss"]] <- "unable, rentrez package missing"
  } else {
    values[["rentrezPkg"]] <- TRUE
  }
  removeModal()
})


observeEvent(input$termButton, {
  validate(
    need(try({
      values[["renInstall"]] == FALSE
    }), "Try again, Press 2.")
  )

  validate(
    need(
      try({
        values[["rentrezPkg"]] == TRUE
      }), values[["renMiss"]]
    )
  )

  showModal(
    modalDialog(
      title = list(
        div("2. Searching string, please wait",
          style = "text-align: center", id = "title0_id"
        ),
        hidden(
          div("Sorry, entrez API seems to be down; Fatal Error",
            style = "text-align: center", id = "title0_id_down"
          )
        )
      ),
      div("this pop-up will close after completion. Press ESC to wait in shiny app",
        style = "text-align: center", id = "subtitle0_id"
      ),
      easyClose = TRUE,
      footer = modalButton(
        list(
          div("Wait in shiny app",
            style = "text-align: center", id = "wait0_id"
          ),
          hidden(
            div("Close",
              class = "title_down",
              style = "text-align: center", id = "wait0_id_down"
            )
          )
        )
      )
    )
  )

  values[["termButtonVal"]] <- 0
  values[["entrezFile"]] <- ""
  values[["rentrezFetch"]] <- 0
  values[["annotated"]] <- TRUE
  values[["button3ab"]] <- NA
  values[["entrez_summary1"]] <- NA
  values[["entrez_titles"]] <- NA
  values[["entrez_selected"]] <- NA
  values[["fetch_listAll"]] <- NA
  values[["authors"]] <- ""

  values[["entrez_search1"]] <- NA
  values[["fetch_list"]] <- NA
  values[["titles_number"]] <- NA

  entrez_summary1 <- list(list(title = numeric(), uid = numeric()))
  entrez_search1 <- list(ids = list())

  term1 <- input$term

  #
  #   search nuccore
  #

  req(input$term)
  entrez_search1 <- tryCatch(rentrez::entrez_search(
      db = "nuccore",
      term = term1,
      retmax = input$maxNum,
      use_history = TRUE
    ),
    error = function(e) {
      "internet or package problem"
    }
  )

  if (entrez_search1[1] != "internet or package problem") {
    if (length(entrez_search1$ids) == 0) {
      values[["searchStatus"]] <- FALSE
      removeModal()
    } else {
      values[["searchStatus"]] <- TRUE
    }
  } else {
    message("internet problem")
    entrez_search1 <- list(ids = list())
    show("title0_id_down")
    show("wait0_id_down")
    hide("subtitle0_id")
    hide("title0_id")
    hide("wait0_id")
    return()
  }

  validate(
    need(try(length(entrez_search1$ids) > 0), "Empty Search")
  )

  values[["entrez_search1"]] <- entrez_search1

  idsvector <- entrez_search1$ids

  lenIds <- length(idsvector)

  if (lenIds > 20) {
    chunkNumber <- lenIds / 20
    idsList <- split(idsvector, sort(rep_len(1:ceiling(chunkNumber), lenIds)))
    entrez_summary1 <- c()
    i <- 1
    while (i <= length(idsList)) {
      Sys.sleep(0.1)
      b <- tryCatch(rentrez::entrez_summary(db = "nuccore", id = idsList[[i]]),
        error = function(e) {
          data.frame(a = "API might be down")
        }
      )
      if (!inherits(b, "data.frame")) {
        entrez_summary1 <- c(entrez_summary1, b)
        i <- i + 1
      } else {
        show("title0_id_down")
        show("wait0_id_down")
        hide("subtitle0_id")
        hide("title0_id")
        hide("wait0_id")
      }
    }
  } else {
    entrez_summary1 <- tryCatch(rentrez::entrez_summary(db = "nuccore", id = idsvector),
      error = function(e) {
        data.frame(a = "API might be down")
      }
    )
    if (inherits(entrez_summary1, "data.frame")) {
      show("title0_id_down")
      show("wait0_id_down")
      hide("subtitle0_id")
      hide("title0_id")
      hide("wait0_id")
    }
  }
  if (!"uid" %in% names(entrez_summary1)) {
    validate(
      need(try(length(entrez_summary1[[1]]$uid) > 0), "wait for summary")
    )
    titles <- character()

    for (i in seq_along(entrez_summary1)) {
      uidTitle <- paste0(i, ". ", entrez_summary1[[i]]$uid, "_", entrez_summary1[[i]]$title)
      titles <- c(titles, uidTitle)
    }
  } else {
    validate(
      need(try(length(entrez_summary1$uid) > 0), "wait for summary")
    )
    titles <- paste0(1, ". ", entrez_summary1$uid, "_", entrez_summary1$title)
  }

  remove(entrez_summary1)
  gc()

  values[["entrez_titles"]] <- titles
  values[["titles_number"]] <- length(titles)

  named_entrez_titles <- seq_along(titles)

  names(named_entrez_titles) <- titles

  values[["named_entrez_titles"]] <- named_entrez_titles

  selection <- grep(term1, titles)

  if (length(selection)) {
    values[["entrez_selected"]] <- selection[1]
  } else {
    values[["entrez_selected"]] <- 1
  }
  values[["termButtonVal"]] <- 1
  removeModal()
})

maxNumReac <- eventReactive(input$termButton, {
  input$maxNum
})

observeEvent(input$button3Download, {
  values[["entrez_selected"]] <- isolate(input$titleSelect)
  showModal(
    modalDialog(
      title = list(
        div("3b. Downloading data, please wait",
          style = "text-align: center", id = "title_id"
        ),
        hidden(
          div("Sorry, entrez API seems to be down; Fatal Error",
            style = "text-align: center", id = "title_id_down"
          )
        ),
        hidden(
          div("Sorry, this doesn't seem a Genome. Unable to proceed!",
            style = "text-align: center", id = "not_a_genome"
          )
        )
      ),
      div("this pop-up will close after completion. Press ESC to wait in shiny app",
        style = "text-align: center", id = "subtitle_id"
      ),
      easyClose = TRUE,
      footer = modalButton(
        list(
          div("Wait in shiny app",
            style = "text-align: center", id = "wait_id"
          ),
          hidden(
            div("Close",
              style = "text-align: center", id = "wait_id_down"
            )
          )
        )
      )
    )
  )

  Sys.sleep(0.4)
  validate(
    need(try(!is.na(values[["entrez_selected"]])), "Empty Search")
  )

  validate(
    need(try(length(values[["entrez_search1"]]$ids) > 0), "Empty Search")
  )

  values[["button3ab"]] <- NA

  values[["names_fetch_list"]] <- NA

  values[["fetch_listAll"]] <- NA
  values[["authors"]] <- ""
  values[["fetch_list"]] <- NA

  values[["geneChrDF"]] <- NA

  mylistAll <- 0
  values[["mylistAllClass"]] <- class(mylistAll)

  mylist <- list()

  values[["entrezFile"]] <- ""
  values[["rentrezFetch"]] <- rentrezFetch <- as.numeric(0)

  idSearch <- values[["entrez_search1"]]$ids[as.numeric(values[["entrez_selected"]])]

  validate(
    need(
      try({
        rentrezFetch <- tryCatch(rentrez::entrez_fetch(
          db = "nuccore",
          id = idSearch,
          rettype = "gbwithparts",
          retmode = "text"
        ), error = function(e) {
          data.frame(a = "API down")
        })
        if (inherits(rentrezFetch, "data.frame")) {
          show("title_id_down")
          show("wait_id_down")
          hide("subtitle_id")
          hide("title_id")
          hide("wait_id")
        }
        class(rentrezFetch) == "character"
      }),
      "downloading"
    )
  )

  values[["entrezFile"]] <- paste0(idSearch, ".gb")
  values[["rentrezFetch"]] <- rentrezFetch

  mylistAll <- tryCatch(genBankReadIF(rentrezFetch), error = function(e) {
    "not a valid format"
  })

  validate(
    need(try({
      if (inherits(mylistAll, "character")) {
        show("not_a_genome")
        show("wait_id_down")
        hide("wait_id")
        hide("subtitle_id")
      }
      class(mylistAll) == "list"
    }), "reading")
  )

  values[["mylistAllClass"]] <- class(mylistAll)

  remove(rentrezFetch)

  values[["authors"]] <- p(paste0(idSearch, ": ", mylistAll$gbdfMain[which(mylistAll$gbdfMain$field == "AUTHORS"), ][1, 2]))

  #
  #   generate d.f. chr. data
  #

  chrDF <- data.frame(chrName = 1, chrSize = mylistAll$source$end)

  mylistAll$source <- NULL

  values[["fetch_listAll"]] <- mylistAll

  chrDF$OTU <- mylistAll$gbdfMain[which(mylistAll$gbdfMain$field == "DEFINITION"), ]$value

  values[["geneChrDF"]] <- chrDF

  genic <- grep("gbdf|source", names(mylistAll), value = TRUE, invert = TRUE)

  if (length(genic) == 0) {
    values[["annotated"]] <- FALSE
    removeModal()
    return()
  } else {
    values[["annotated"]] <- TRUE
  }

  mylist <- mylistAll[which(names(mylistAll) %in% genic)]

  values[["names_fetch_list"]] <- names(mylist)
  values[["fetch_list_selected"]] <- grep("CDS", names(mylist), invert = TRUE)
  values[["fetch_list"]] <- mylist

  validate(
    need(try(!invalid(values[["names_fetch_list"]])), "downloading 804")
  )

  values[["button3ab"]] <- as.numeric(input$button3Download)
  gc()
  removeModal()
})

observeEvent(input$makeDFsButton, {
  showModal(modalDialog(
    title = "4. Making data.frames, please wait",
    "this pop-up will close after completion, Press ESC to wait in shiny app",
    easyClose = TRUE,
    footer = modalButton("Wait in shiny app")
  ))

  values[["geneMarkDF"]] <- values[["geneMarkDFOrig"]] <- NULL

  mylist <- values[["fetch_list"]]

  chrDF <- values[["geneChrDF"]]

  mylist <- mylist[as.numeric(input$fetchSelect)]

  mylistChrDF <- dplyr::bind_rows(mylist, .id = "feature")
  # add necessary columns
  mylistChrDF$markPos <- pmin(as.numeric(mylistChrDF$begin), as.numeric(mylistChrDF$end))
  mylistChrDF$markSize <- abs(as.numeric(mylistChrDF$end) - as.numeric(mylistChrDF$begin))

  #
  #   MARK NAMES
  mylistChrDF$markName <- mylistChrDF$locus_tag

  # careful
  # mylistChrDF[which(is.na(mylistChrDF$markName) ),]$markName<-
  # sub("([[:alpha:] ]+);.*","\\1", mylistChrDF[which(is.na(mylistChrDF$markName) ),]$note )

  # orientation of arrows
  mylistChrDF$style <- ifelse(mylistChrDF$isComplement, "downArrow", "upArrow")

  mylistChrDF$chrName <- unique(chrDF$chrName)
  # select main columns for data.frame of marks' positions
  columnList <- c(
    "chrName", "markName", "markPos", "markSize",
    "feature", "isJoin", "style", "pseudo", "gene",
    "regulatory_class"
  )

  columnMarkList <- intersect(columnList, colnames(mylistChrDF))

  # remove some columns
  marksDfChr <- mylistChrDF[, columnMarkList]

  marksDfChr$OTU <- chrDF$OTU
  # add mandatory column

  rownames(marksDfChr) <- 1:nrow(marksDfChr)

  values[["pseudo"]] <- 0

  if ("pseudo" %in% colnames(marksDfChr)) {
    values[["pseudo"]] <- nrow(marksDfChr[which(marksDfChr$pseudo == TRUE), ])
  }

  values[["geneMarkDF"]] <- values[["geneMarkDFOrig"]] <- marksDfChr

  markStyle <- makedfMarkColorMycolors(
    unique(marksDfChr$markName)
    # , c("black","forestgreen","cornflowerblue")
    , if (length(values[["mycolors2"]]) == 0) {
      ""
    } else {
      values[["mycolors2"]]
    }
  )

  # arrows
  markStyle$style <- marksDfChr$style[match(markStyle$markName, marksDfChr$markName)]
  markStyle$protruding <- marksDfChr$protruding[match(markStyle$markName, marksDfChr$markName)]

  values[["markStyleDF"]] <- markStyle
  removeModal()
  gc()
})


observeEvent(input$modifyMarksButton, {
  showModal(modalDialog(
    title = "5. Modifying marks, please wait",
    "this pop-up will close after completion, Press ESC to wait in shiny app",
    easyClose = TRUE,
    footer = modalButton("Wait in shiny app")
  ))

  values[["geneMarkDF"]] <- NA
  values[["markStyleDF"]] <- NA
  values[["geneMarkDFOrigCopy"]] <- NA
  values[["geneMarkDFOrigCopy"]] <- values[["geneMarkDFOrig"]]

  #
  #   modify mark names in columns
  #

  if (input$pseudo == "onlyPseudo" &
    "pseudo" %in% colnames(values[["geneMarkDFOrigCopy"]])) {
    values[["geneMarkDFOrigCopy"]] <- values[["geneMarkDFOrigCopy"]][which(values[["geneMarkDFOrigCopy"]]$pseudo == TRUE), ]
  } else if (input$pseudo == "removePseudo" &
    "pseudo" %in% colnames(values[["geneMarkDFOrigCopy"]])) {
    values[["geneMarkDFOrigCopy"]] <- values[["geneMarkDFOrigCopy"]][which(is.na(values[["geneMarkDFOrigCopy"]]$pseudo)), ]
  }

  # Replace codes with genes, and replace NAs in markNames (locus_tag)
  if (input$useGeneNames) {
    tryCatch(values[["geneMarkDFOrigCopy"]][which(!is.na(values[["geneMarkDFOrigCopy"]]$gene)), ]$markName <-
      values[["geneMarkDFOrigCopy"]][which(!is.na(values[["geneMarkDFOrigCopy"]]$gene)), ]$gene, error = function(e) {
      ""
    })
  }

  if (input$useRCNames) {
    tryCatch(values[["geneMarkDFOrigCopy"]][which(!is.na(values[["geneMarkDFOrigCopy"]]$regulatory_class)), ]$markName <-
      values[["geneMarkDFOrigCopy"]][which(!is.na(values[["geneMarkDFOrigCopy"]]$regulatory_class)), ]$regulatory_class, error = function(e) {
      ""
    })
  }
  if (input$makeUnique) {
    values[["geneMarkDFOrigCopy"]]$markNameOrig <- values[["geneMarkDFOrigCopy"]]$markName
    values[["geneMarkDFOrigCopy"]]$markName <- make.uniqueIF(values[["geneMarkDFOrigCopy"]]$markNameOrig)
  }

  marksDfChrCols <- namesToColumns(values[["geneMarkDFOrigCopy"]],
    values[["geneChrDF"]],
    markType = input$markType # c("downArrow"),
    , amountofSpaces = input$amountofSpaces # 10
    , colNumber = input$colNumber # 2,
    , protrudingInt = input$protrudingInt # 1.3
    , protruding = input$protruding,
    circularPlot = input$circularPlot,
    rotation = input$rotation
  )

  if (input$addSTARTPos) {
    marksDfChrCols <- plyr::rbind.fill(
      marksDfChrCols,
      data.frame(
        markName = paste0(
          paste0(rep(" ", input$colNumber * input$amountofSpaces),
            collapse = ""
          ),
          "START"
        ),
        markPos = 1,
        markSize = NA,
        style = "square",
        OTU = unique(marksDfChrCols$OTU),
        chrName = unique(marksDfChrCols$chrName)
      )
    )
  }

  if (input$nucMarkStyle == "cM") {
    tryCatch(marksDfChrCols[which(marksDfChrCols$style %in% "downArrow"), ]$style <- "cMLeft", error = function(e) {})
    tryCatch(marksDfChrCols[which(marksDfChrCols$style %in% "upArrow"), ]$style <- "cM", error = function(e) {})
  }

  # create mark general data data.frame
  markStyle <- makedfMarkColorMycolors(
    unique(marksDfChrCols$markName),
    if (length(values[["mycolors2"]]) == 0) {
      ""
    } else {
      values[["mycolors2"]]
    }
  )

  # arrows
  markStyle$style <- marksDfChrCols$style[match(markStyle$markName, marksDfChrCols$markName)]
  markStyle$protruding <- marksDfChrCols$protruding[match(markStyle$markName, marksDfChrCols$markName)]


  if (input$colorFeature) {
    markStyle$feature <- as.character(marksDfChrCols$feature[match(markStyle$markName, marksDfChrCols$markName)])
    markStyleFeature <- makedfMarkColorMycolors(
      unique(marksDfChrCols$feature),
      if (length(values[["mycolors2"]]) == 0) {
        ""
      } else {
        values[["mycolors2"]]
      }
    )
    markStyle$markColor <- markStyleFeature$markColor[match(markStyle$feature, markStyleFeature$markName)]
  }

  if (input$mirror) {
    marksDfChrCols$markPos <- as.numeric(values[["geneChrDF"]]$chrSize) - as.numeric(marksDfChrCols$markPos) - as.numeric(marksDfChrCols$markSize)

    tryCatch(markStyle$style[which(markStyle$style == "downArrow")] <- "upArrow1", error = function(e) {
      ""
    })
    tryCatch(markStyle$style[which(markStyle$style == "cMLeft")] <- "cM1", error = function(e) {
      ""
    })

    tryCatch(markStyle$style[which(markStyle$style == "upArrow")] <- "downArrow", error = function(e) {
      ""
    })
    tryCatch(markStyle$style[which(markStyle$style == "cM")] <- "cMLeft", error = function(e) {
      ""
    })

    tryCatch(markStyle$style[which(markStyle$style == "upArrow1")] <- "upArrow", error = function(e) {
      ""
    })
    tryCatch(markStyle$style[which(markStyle$style == "cM1")] <- "cM", error = function(e) {
      ""
    })
  }

  values[["geneMarkDF"]] <- marksDfChrCols
  values[["markStyleDF"]] <- markStyle

  removeModal()
  gc()
})

observeEvent(input$loadDFbutton, {
  showModal(modalDialog(
    title = "6. Loading data.frames, please wait",
    "this pop-up will close after completion (See data.frames page), Press ESC to wait in shiny app",
    easyClose = TRUE,
    footer = modalButton("Wait in shiny app")
  ))

  values[["df1"]] <- data.frame()
  values[["df1Mark"]] <- data.frame()
  values[["df1MStyle"]] <- data.frame()
  Sys.sleep(3)

  validate(need(try(values[["geneChrDF"]]), "not ready chr. d.f."))
  validate(need(try(inherits(values[["geneChrDF"]], "data.frame")), "still not ready chr. d.f."))
  values[["df1"]] <- values[["geneChrDF"]]

  validate(need(try(values[["geneMarkDF"]]), "not ready mark d.f."))
  validate(need(try(inherits(values[["geneMarkDF"]], "data.frame")), "still not ready mark d.f."))
  values[["df1Mark"]] <- values[["geneMarkDF"]]
  values[["geneMarkDF"]] <- NULL

  validate(need(try(values[["markStyleDF"]]), "not ready mark style d.f."))
  validate(need(try(inherits(values[["markStyleDF"]], "data.frame")), "still not ready mark style d.f."))
  values[["df1MStyle"]] <- values[["markStyleDF"]]
  values[["markStyleDF"]] <- NULL
  removeModal()
})



observeEvent(values[["pseudo"]], {
  updateRadioButtons(session, "pseudo", label = paste0("Show: ", values[["pseudo"]], " found"))
})
