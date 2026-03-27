## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev = "ragg_png",
  fig.width = 6, 
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

## ----plain, fig.cap="", fig.alt="Example of a plain ggplot2"------------------
ggplot(df, aes(trt, outcome)) +
  geom_col(aes(fill=trt),colour='black') +
  theme_bw() +
  labs(title = "Plain ggplot2")

## ----ggpattern, fig.cap="", fig.alt="Example using ggpattern to make patterned geoms"----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(aes(fill=trt),colour='black', pattern = 'circle') +
  theme_bw() +
  labs(title = "ggpattern") +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----density, fig.cap="", fig.alt="Example of increasing density of striping"----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill=trt, pattern_shape=trt),
    colour          = 'black', 
    pattern         = 'pch',
    pattern_density = 0.5
  ) +
  theme_bw() +
  labs(title = "Fixed density of 0.5 (50% of the fill area)") + 
  theme(legend.key.size = unit(1.5, 'cm'))

## ----density2, fig.cap="", fig.alt="Example from using the density aesthetic as a mapped aesthetic"----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_density = trt),
    colour          = 'black', 
    pattern         = 'circle'
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Density") + 
  theme(legend.key.size = unit(1.5, 'cm'))

## ----density3, fig.cap="", fig.alt="Example from using the density aesthetic as a mapped aesthetic with manual scale."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_density = trt),
    colour          = 'black', 
    pattern         = 'regular_polygon',
    pattern_shape   = 'convex6',
    pattern_fill    = 'white',
    pattern_colour  = NA,
    pattern_grid    = 'hex'
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Density") + 
  theme(legend.key.size = unit(1.5, 'cm')) + 
  scale_pattern_density_manual(values = c(a = 0.7, b=0.8, c=0.9))

## ----spacing, fig.cap="", fig.alt="Example from using the spacing aesthetic as a mapped aesthetic."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_spacing = trt, pattern_shape = trt),
    pattern_density = 0.7,
    colour          = 'black', 
    pattern         = 'regular_polygon',
    pattern_scale   = 0.2
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Spacing") + 
  scale_pattern_shape_manual(values = c(a = "circle", b = "convex3", c = "star5")) +
  theme(legend.key.size = unit(1.5, 'cm'))

## ----fill, fig.cap="", fig.alt="Example from using the fill aesthetic as a mapped aesthetic."----
ggplot(df, aes(trt, outcome)) +
  geom_col_pattern(
    aes(fill = trt, pattern_fill = trt),
    colour          = 'black', 
    pattern         = 'rose',
    pattern_frequency = 5/2,
    pattern_density = 0.7,
    pattern_spacing = 0.15
  ) +
  theme_bw() +
  labs(title = "Aesthetic Mapping of 'trt' to Pattern Fill") + 
  scale_pattern_fill_viridis_d() + 
  theme(legend.key.size = unit(1.5, 'cm'))

