## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = F-----------------------------------------------------------------
#  install.packages("qqconf")

## ----message = F--------------------------------------------------------------
require(qqconf)

## -----------------------------------------------------------------------------
set.seed(0)
sample <- rnorm(n = 100, mean = 0, sd = 1)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1)
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1),
  ylim = c(-4, 4)
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1),
  ylim = c(-4, 4),
  xlab = "More Informative Title"
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1),
  ylim = c(-4, 4),
  points_params = list(cex = .5) # makes points smaller
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1),
  ylim = c(-4, 4),
  line_params = list(col="red") # makes expectation line red
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1),
  ylim = c(-4, 4),
  line_params = list(type="n") # removes expectation line
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
qq_conf_plot(
  obs = sample, 
  distribution = qnorm,
  dparams = list(mean = 0, sd = 1),
  ylim = c(-4, 4),
  polygon_params = list(col = 'powderblue', border = NA) # change shading and keep no border
)

## -----------------------------------------------------------------------------
sample <- rbeta(n = 10000, shape1 = 1.0, shape2 = 1.05)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
pp_conf_plot(
  obs = sample, 
  distribution = punif,
  points_params = list(cex=.1)
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
pp_conf_plot(
  obs = sample, 
  distribution = punif,
  points_params = list(cex=.1),
  difference = TRUE, # Make y-axis differenced
  ylim = c(-.0225, .0225)
)

## -----------------------------------------------------------------------------
mix <- distr::UnivarMixingDistribution(
  distr::Beta(shape1 = .25, shape2 = 1), 
  distr::Unif(),
  mixCoeff=c(.01, .99)
)

sampler <- distr::r(mix)

sample <- sampler(10000)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
pp_conf_plot(
  obs = sample, 
  distribution = punif,
  points_params = list(cex=.1)
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
pp_conf_plot(
  obs = sample, 
  distribution = punif,
  points_params = list(cex=.1),
  difference = TRUE
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
pp_conf_plot(
  obs = sample, 
  distribution = punif,
  points_params = list(cex=.1),
  log10 = TRUE
)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
bounds <- get_bounds_two_sided(alpha = .05, n = 100)

## ----fig.align = "center", fig.width = 7, fig.height = 5----------------------
bounds <- get_bounds_one_sided(alpha = .05, n = 100)

## ----session-info-------------------------------------------------------------
sessionInfo()

