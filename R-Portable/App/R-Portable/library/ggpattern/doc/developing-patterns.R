## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
dev = "ragg_png",
comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("ggplot2")
library("ggpattern")

## ----echo = FALSE-------------------------------------------------------------
x <- read.csv(textConnection("
x,  y,  id
0,  0,   1
1,  0,   1
1,  1,   1
0,  1,   1
0,  0,   2
2,  0,   2
2,  1,   2
0,  1,   2"))

knitr::kable(x, caption = "example data in 'polygon_df' format")

## ----eval = FALSE-------------------------------------------------------------
# options(ggpattern_array_funcs    = list(your_pattern_name = your_pattern_function))
# options(ggpattern_geometry_funcs = list(your_pattern_name = your_pattern_function))

## -----------------------------------------------------------------------------
tiling3_pattern <- function(params, boundary_df, aspect_ratio, legend = FALSE) {
    args <- as.list(params)
    args <- args[grep("^pattern_", names(args))]

    # hexagonal tiling using "regular_polygon" pattern
    args$pattern <- "polygon_tiling"

    # three-color tiling using `fill`, `pattern_fill` and their "average"
    avg_col <- gridpattern::mean_col(params$fill, params$pattern_fill)
    args$pattern_fill <- c(params$fill, avg_col, args$pattern_fill)

    args$x <- boundary_df$x
    args$y <- boundary_df$y
    args$id <- boundary_df$id
    args$prefix <- ""

    do.call(gridpattern::patternGrob, args)
}


## -----------------------------------------------------------------------------
options(ggpattern_geometry_funcs = list(tiling3 = tiling3_pattern))

## ----tiling3_pattern, fig.alt = "ggplot2 plot using a custom polygon tiling pattern fill."----
df <- data.frame(trt = c("a", "b", "c"), outcome = c(2.3, 1.9, 3.2))
ggplot(df, aes(trt, outcome)) +
    geom_col_pattern(aes(fill = trt, pattern_type = trt),
                     pattern = 'tiling3',
                     pattern_angle = 45,
                     pattern_spacing = 0.15) +
    scale_pattern_type_manual(values = c("hexagonal", "tetrakis_square", "rhombille")) +
    theme(legend.key.size = unit(1.2, 'cm'))

