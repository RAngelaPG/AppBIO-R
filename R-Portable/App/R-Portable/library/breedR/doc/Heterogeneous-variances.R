## ----setup, echo = FALSE, include=FALSE---------------------------------------
library(breedR)
library(ggplot2)
theme_set(theme_bw())

## ----sim-data-----------------------------------------------------------------
set.seed(123)

n <- 1e3   # n obs
sigma2 <- 4  # true residual variance (for a weight of 1)
w = runif(n, min = .5, max = 2)  # vector of weights

dat <- 
  transform(
    data.frame(
      e = rnorm(n, sd = sqrt(sigma2))
    ),
    y = 10 + e/sqrt(w)  # simulated phenotype
  )

## ----wheights-example---------------------------------------------------------
res <- remlf90(
  y ~ 1,
  data = dat,
  weights = w  # specification of weights
)

## ----weights-results----------------------------------------------------------
summary(res)

ggplot(transform(dat, est_e = residuals(res)), aes(e, est_e)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "darkgray")

