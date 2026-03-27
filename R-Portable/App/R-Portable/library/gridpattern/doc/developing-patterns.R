## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----setup--------------------------------------------------------------------
library("gridpattern")

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
create_pattern_polygon <- function(params, boundary_df, aspect_ratio, legend = FALSE) {
    x <- boundary_df$x
    y <- boundary_df$y
    id <- boundary_df$id

    alpha <- ifelse(is.na(params$pattern_alpha), 1, params$pattern_alpha)
    gp <- grid::gpar(alpha = alpha, 
                     col = params$pattern_colour,
                     fill = params$pattern_fill,
                     lty = params$pattern_linetype,
                     lwd = params$pattern_linewidth)

    grid::polygonGrob(x = x, y = y, id = id, default.units = "npc", gp = gp)
}

## -----------------------------------------------------------------------------
options(ggpattern_geometry_funcs = list(polygon = create_pattern_polygon))

## ----fig.alt = "Three rectangles each filled with a single color"-------------
grid.pattern("polygon", fill = "red", size = 4, linetype = "dashed", 
             x = c(0.05, 0.05, 0.305, 0.305), y = c(0.05, 0.305, 0.305, 0.05))
grid.pattern("polygon", fill = "green", alpha = 0.2, 
             x = c(0.35, 0.35, 0.65, 0.65), y = c(0.35, 0.65, 0.65, 0.35))
grid.pattern("polygon", fill = "blue", colour = "grey",
             x = c(0.7, 0.7, 1.0, 1.0), y = c(0.7, 1.0, 1.0, 0.7))

## -----------------------------------------------------------------------------
create_pattern_complex <- function(params, boundary_df, aspect_ratio, legend = FALSE) {
    args <- as.list(params)
    args <- args[grep("^pattern_", names(args))]
    args$x <- boundary_df$x
    args$y <- boundary_df$y
    args$id <- boundary_df$id
    args$prefix <- ""

    args_stripe <- args
    args_stripe$pattern <- "stripe"
    args_stripe$pattern_density <- 0.5 * args$pattern_density
    args_stripe$pattern_spacing <- 2 * args$pattern_spacing
    grob_stripe <- do.call(gridpattern::patternGrob, args_stripe)

    args_circle <- args
    args_circle$pattern <- "regular_polygon"
    args_circle$pattern_shape <- c("circle", "null")
    args_circle$pattern_yoffset <- args$pattern_spacing + args$pattern_yoffset
    args_circle$pattern_type = "horizontal"
    grob_circle <- do.call(gridpattern::patternGrob, args_circle)

    args_gradient <- args
    args_gradient$pattern <- "gradient"
    args_gradient$pattern_fill <- "#00000070"
    args_gradient$pattern_fill2 <- "#FFFFFF70"
    args_gradient$pattern_orientation <- "vertical"
    grob_gradient <- do.call(gridpattern::patternGrob, args_gradient)

    grid::grobTree(grob_stripe, grob_circle, grob_gradient)
}

## -----------------------------------------------------------------------------
options(ggpattern_geometry_funcs = list(complex = create_pattern_complex))

## ----fig.alt = "Three rectangles each with a 'stripe', 'circle', and 'gradient' pattern"----
grid.pattern("complex", fill = "red", angle = 45, spacing = 0.05, density = 0.3,
             x = c(0.0, 0.0, 0.3, 0.3), y = c(0.0, 0.3, 0.3, 0.0))
grid.pattern("complex", fill = "green", angle = 45, spacing = 0.2, density = 0.2,
             x = c(0.35, 0.35, 0.65, 0.65), y = c(0.35, 0.65, 0.65, 0.35))
grid.pattern("complex", fill = "blue", angle = 45, spacing = 0.1, density = 0.3,
             x = c(0.7, 0.7, 1.0, 1.0), y = c(0.7, 1.0, 1.0, 0.7))
grid::grid.polygon(x = c(0.0, 0.0, 0.3, 0.3, 0.35, 0.35, 0.65, 0.65, 0.7, 0.7, 1.0, 1.0),
                   y = c(0.0, 0.3, 0.3, 0.0, 0.35, 0.65, 0.65, 0.35, 0.7, 1.0, 1.0, 0.7), 
                   id = rep(1:3, each = 4), 
                   gp = grid::gpar(col = "black", fill = NA, lwd=4))

## -----------------------------------------------------------------------------
create_pattern_simple <- function(width, height, params, legend) {
  
  # Ensure the selected pattern is sane.
  choice <- params$pattern_type
  if (is.null(choice) || is.na(choice) || !is.character(choice)) {
    choice <- 'a'
  }
  
  # Choose the values with which to fill the array
  values <- switch(
    choice,
    a = rep(c(0, 1, 0, 1, 1, 0, 0, 1, 1, 1), each = 3),
    b = rep(c(1, 0, 0, 1, 0.5, 0.5, 1, 1, 0, 0, 0, 0, 0, 0.5), each = 7),
    c = rep(seq(0, 1, 0.05), each = 7),
    rep(c(0, 1, 0, 1, 1, 0, 0, 1, 1, 1), each = 3)
  )
  # Create an RGBA array of the requested dimensions
  simple_array <- array(values, dim = c(height, width, 4))

  simple_array
}

## -----------------------------------------------------------------------------
options(ggpattern_array_funcs = list(simple = create_pattern_simple))

## ----fig.alt = "Three rectangles each filled with an 'array' pattern"---------
grid::grid.polygon(x = c(0, 0, 1, 1), y = c(0, 1, 1, 0), 
                   gp = grid::gpar(col=NA, fill="grey"))
grid.pattern("simple", type = "a",
             x = c(0.0, 0.0, 0.3, 0.3), y = c(0.0, 0.3, 0.3, 0.0))
grid.pattern("simple", type = "b",
             x = c(0.35, 0.35, 0.65, 0.65), y = c(0.35, 0.65, 0.65, 0.35))
grid.pattern("simple", type = "c",
             x = c(0.7, 0.7, 1.0, 1.0), y = c(0.7, 1.0, 1.0, 0.7))

