## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----setup--------------------------------------------------------------------
#  library(shinylogs)

## -----------------------------------------------------------------------------
#  # server
#  track_usage(storage_mode = store_json(path = "logs/"))

## -----------------------------------------------------------------------------
#  # UI
#  use_tracking()

## ----ex-session, eval=TRUE, echo=FALSE----------------------------------------
data.frame(stringsAsFactors=FALSE,
                   app = c("iris-cluster", "iris-cluster", "iris-cluster", "iris-cluster", "iris-cluster", "iris-cluster"),
                  user = c("dreamRs", "dreamRs", "dreamRs", "dreamRs", "dreamRs", "dreamRs"),
      server_connected = c("2019-06-19 14:07:09", "2019-06-19 14:56:32",
                           "2019-06-19 14:57:45", "2019-06-19 15:03:00",
                           "2019-06-19 15:08:53", "2019-08-04 10:58:47"),
             sessionid = c("9815194cfaa6317fbea68ae9537d63d1",
                           "0feeacf3201f5cca088059cec2b0a710",
                           "8b12c138f159cc5805e136338dddac59",
                           "f3950a1588b78d45c53c756a4136df67", "254cafb3d5866384f13022d066546cae",
                           "8a431f4b90b3b1926047dd2539be0793"),
   server_disconnected = c("2019-06-19 14:07:17", "2019-06-19 14:57:39",
                           "2019-06-19 15:01:36", "2019-06-19 15:08:13",
                           "2019-06-19 15:09:38", "2019-08-04 10:58:49"),
            user_agent = c("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36",
                           "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36",
                           "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36",
                           "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36",
                           "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/74.0.3729.169 Safari/537.36", NA),
            screen_res = c("1920x1080", "1920x1080", "1920x1080", "1920x1080",
                           "1920x1080", NA),
           browser_res = c("1574x724", "1920x937", "1920x937", "1920x937",
                           "1920x937", NA),
           pixel_ratio = c(1L, 1L, 1L, 1L, 1L, NA),
     browser_connected = c("2019-06-19 14:07:09", "2019-06-19 14:56:32",
                           "2019-06-19 14:57:45", "2019-06-19 15:03:00",
                           "2019-06-19 15:08:53", NA)
)

## ----ex-inputs, eval=TRUE, echo=FALSE-----------------------------------------
data.frame(stringsAsFactors=FALSE,
   sessionid = c("9815194cfaa6317fbea68ae9537d63d1",
                 "9815194cfaa6317fbea68ae9537d63d1",
                 "9815194cfaa6317fbea68ae9537d63d1", "9815194cfaa6317fbea68ae9537d63d1",
                 "9815194cfaa6317fbea68ae9537d63d1", "9815194cfaa6317fbea68ae9537d63d1"),
        name = c("xcol", "clusters", "clusters", "clusters", "clusters",
                 "clusters"),
   timestamp = c("2019-06-19 14:07:11", "2019-06-19 14:07:15",
                 "2019-06-19 14:07:15", "2019-06-19 14:07:15",
                 "2019-06-19 14:07:14", "2019-06-19 14:07:13"),
       value = c("Sepal.Width", 6, 4, 5, 7, 4),
        type = c(NA, "shiny.number", "shiny.number", "shiny.number",
                 "shiny.number", "shiny.number"),
     binding = c("shiny.selectInput", "shiny.numberInput", "shiny.numberInput",
                 "shiny.numberInput", "shiny.numberInput",
                 "shiny.numberInput")
)

## ----ex-outputs, eval=TRUE, echo=FALSE----------------------------------------
data.frame(stringsAsFactors=FALSE,
   sessionid = c("9815194cfaa6317fbea68ae9537d63d1",
                 "9815194cfaa6317fbea68ae9537d63d1",
                 "9815194cfaa6317fbea68ae9537d63d1", "9815194cfaa6317fbea68ae9537d63d1",
                 "9815194cfaa6317fbea68ae9537d63d1", "9815194cfaa6317fbea68ae9537d63d1"),
        name = c("plot1", "plot1", "plot1", "plot1", "plot1", "plot1"),
   timestamp = c("2019-06-19 14:07:15", "2019-06-19 14:07:15",
                 "2019-06-19 14:07:14", "2019-06-19 14:07:15",
                 "2019-06-19 14:07:11", "2019-06-19 14:07:10"),
     binding = c("shiny.imageOutput", "shiny.imageOutput", "shiny.imageOutput",
                 "shiny.imageOutput", "shiny.imageOutput",
                 "shiny.imageOutput")
)

## ----ex-errors, eval=TRUE, echo=FALSE-----------------------------------------
data.frame(stringsAsFactors=FALSE,
   sessionid = c("f8f50a3743023aae7d0d6350a2fd6841"),
        name = c("plot1"),
   timestamp = c("2019-06-19 14:07:18"),
       error = c("NA/NaN/Inf in foreign function call (arg 1)")
)

