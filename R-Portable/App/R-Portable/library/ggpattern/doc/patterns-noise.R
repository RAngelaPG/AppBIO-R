## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "ragg_png",
  fig.width = 8,
  fig.height = 6
)

# Can't use {ragg} to save gradient example on CRAN due to UBSAN issue
if (!identical(Sys.getenv("IN_PKGDOWN"), "true")) {
  knitr::opts_chunk$set(
    dev = "png",
    dev.args = list(type = "cairo"),
    eval = all(capabilities(c("png", "cairo")))
  )
}

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages({
  library("ggplot2")
  library("ggpattern")
  require("magick", quietly = TRUE)
})

## -----------------------------------------------------------------------------
df <- data.frame(
  trt     = c("a", "b", "c"),
  outcome = c(2.3, 1.9, 3.2)
)

## ----plasma, fig.cap="", fig.alt="Example of 'plasma' pattern"----------------
if (require("magick")) {

ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(
      fill         = trt,
      pattern_fill = trt
    ),
    fill          = NA,
    pattern       = 'plasma',
    pattern_alpha = 1,
    pattern_scale = 2,
    colour        = 'black'
  ) +
  theme_bw(15) +
  labs(
    title    = "ggpattern::geom_col_pattern()",
    subtitle = "pattern='plasma'"
  ) +
  theme(legend.key.size = unit(1.5, 'cm'))

}

## ----plasma2, fig.cap="", fig.alt="Example of 'plasma' pattern with 'pattern_alpha = 0.7'"----
if (require("magick")) {

ggplot(mtcars) +
  geom_density_pattern(
    aes(
      x             = mpg,
      pattern_fill  = as.factor(cyl)
    ),
    pattern      = 'plasma',
    pattern_alpha = 0.7
  ) +
  theme_bw(15) +
  theme(legend.position = 'none') +
  labs(
    title    = "ggpattern::geom_density_pattern()",
    subtitle = "pattern='plasma'"
  )

}

## ----gradient, fig.cap="", fig.alt="Example of 'gradient' pattern fade to white"----
if (require("magick")) {

ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(
      pattern_fill        = trt,
      pattern_orientation = trt
    ),
    pattern       = 'gradient',
    pattern_fill2 = 'white',
    colour        = 'black'
  ) +
  theme_bw(15) +
  labs(
    title    = "ggpattern::geom_col_pattern()",
    subtitle = "pattern = 'gradient'"
  ) +
  theme(
    legend.key.size = unit(1.5, 'cm')
  )

}

## ----gradient2, fig.cap="", fig.alt="Example of 'gradient' pattern fade to dark blue"----
if (require("magick")) {

ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(
      pattern_fill        = trt,
      pattern_orientation = trt
    ),
    pattern       = 'gradient',
    pattern_fill2 = '#445566',
    colour        = 'black'
  ) +
  theme_bw(15) +
  labs(
    title    = "ggpattern::geom_col_pattern()",
    subtitle = "pattern = 'gradient'"
  ) +
  theme(
    legend.key.size = unit(1.5, 'cm')
  )

}

## ----gradient3, fig.cap="", fig.alt="Example of 'gradient' pattern fade to transparent"----
if (require("magick")) {

ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(
      pattern_fill        = trt,
      pattern_orientation = trt
    ),
    pattern       = 'gradient',
    pattern_fill2 = NA,
    fill          = NA,
    colour        = 'black'
  ) +
  theme_bw(15) +
  labs(
    title    = "ggpattern::geom_col_pattern()",
    subtitle = "pattern = 'gradient'"
  ) +
  theme(legend.key.size = unit(1.5, 'cm'))

}

## ----gradient4, fig.cap="", fig.alt="Example of 'gradient' pattern with a non-rectangular geom"----
if (require("magick")) {

ggplot(mtcars) +
  geom_density_pattern(
    aes(
      x = mpg,
      pattern_fill = as.factor(cyl),
      pattern_orientation = as.factor(cyl)
    ),
    pattern       = 'gradient',
    pattern_fill2 = NA,
    fill          = NA
  ) +
  theme_bw(15) +
  labs(
    title    = "ggpattern::geom_density_pattern()",
    subtitle = "pattern = 'gradient'"
  ) +
  theme(legend.key.size = unit(1.5, 'cm'))

}

## ----ambient, fig.cap="", fig.alt="Example of 'ambient' pattern"--------------
if (require("ambient")) {
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(pattern_fill = trt),
    pattern       = 'ambient',
    pattern_fill2 = 'white',
    colour        = NA,
    fill          = NA
  ) +
  theme_bw(15) +
  labs(
    title    = "ggpattern::geom_density_pattern()",
    subtitle = "pattern = 'ambient'"
  ) +
  theme(legend.position = 'none')
}

## ----ambient2, fig.cap="", fig.alt="Another example of 'ambient' pattern"-----
if (require("ambient")) {
ggplot(mtcars) +
  geom_density_pattern(
    aes(
      x = mpg,
      pattern_fill  = as.factor(cyl),
      pattern_fill2 = as.factor(cyl)
    ),
    pattern = 'ambient'
  ) +
  theme_bw(15) +
  theme(legend.key.size = unit(2, 'cm')) +
  scale_pattern_fill_brewer (palette = 'Accent', direction =  1) +
  scale_pattern_fill2_brewer(palette = 'Dark2' , direction =  1) +
  labs(
    title    = "ggpattern::geom_density_pattern()",
    subtitle = "pattern = 'ambient'"
  )
}

