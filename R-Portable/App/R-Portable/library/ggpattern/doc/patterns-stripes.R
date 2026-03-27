## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "ragg_png",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpattern)
})

## -----------------------------------------------------------------------------
df <- data.frame(trt = c("a", "b", "c"), outcome = c(2.3, 1.9, 3.2))
df

## ----vanilla, fig.cap="", fig.alt="A plot using vanilla ggplot2"--------------
ggplot(df, aes(trt, outcome)) +
  geom_col(aes(fill=trt),colour='black') +
  theme_bw() +
  labs(title = "Plain ggplot2")

## ----ggpattern, fig.cap="", fig.alt="A plot using ggpattern with a patterned geom"----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(aes(fill=trt),colour='black',
                   pattern='stripe') +
  theme_bw() +
  labs(title = "ggpattern") +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----mixing, fig.cap="", fig.alt="A plot with both a 'stripe' pattern and a 'wave' pattern."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(aes(fill=trt, pattern=trt, pattern_type=trt),
                   colour='black', pattern_key_scale_factor=0.5,
                   pattern_spacing=0.7, pattern_frequency=0.7, pattern_units='cm') +
  theme_bw() +
  labs(title = "Use 'stripe' and 'wave' patterns") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  scale_pattern_manual(values=c('stripe', 'wave', 'wave')) +
  scale_pattern_type_manual(values=c(NA, 'triangle', 'sine'))

## ----density, fig.cap="", fig.alt="A plot with the density of striping increased to 50% of fill area"----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill=trt, pattern=trt),
    colour          = 'black',
    pattern_density = 0.5
  ) +
  theme_bw() +
  labs(title = "Fixed density of 0.5 (50% of the fill area)") +
  scale_pattern_manual(values=c('stripe', 'crosshatch', 'weave')) +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----weave, fig.cap="", fig.alt="A plot with 'plain', 'twill', 'satin' weave patterns"----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(pattern_fill2=trt, pattern_type=trt),
    pattern = 'weave',
    colour          = 'black',
    pattern_density = 1.0,
    pattern_fill = 'grey',
    pattern_key_scale_factor = 0.5,
  ) +
  theme_bw() +
  labs(title = "Some 'weave' types") +
  scale_pattern_type_manual(values=c('plain', 'twill', 'satin')) +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----density2, fig.cap="", fig.alt="Example plot from using the density aesthetic as a mapped aesthetic."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_density = trt),
    colour          = 'black',
    pattern         = 'stripe'
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Density") +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----density3, fig.cap="", fig.alt="Example plot from using the density aesthetic as a mapped aesthetic with manual scale."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_density = trt),
    colour          = 'black',
    pattern         = 'stripe'
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Density") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  scale_pattern_density_manual(values = c(a = 0.1, b=0.3, c=0.5))

## ----spacing, fig.cap="", fig.alt="Example plot from using the spacing aesthetic as a mapped aesthetic."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_spacing = trt),
    colour          = 'black',
    pattern         = 'stripe'
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Spacing") +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----fill, fig.cap="", fig.alt="Example plot from using the fill aesthetic as a mapped aesthetic."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_fill = trt),
    colour          = 'black',
    pattern         = 'stripe'
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Pattern Fill") +
  scale_pattern_fill_viridis_d() +
  theme(legend.key.size = unit(1.5, 'cm'))

