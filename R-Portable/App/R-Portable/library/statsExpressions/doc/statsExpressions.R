## -----------------------------------------------------------------------------
source("setup.R")

## -----------------------------------------------------------------------------
citation("statsExpressions")

## -----------------------------------------------------------------------------
mtcars %>% oneway_anova(cyl, wt, type = "nonparametric")

mtcars %>% oneway_anova(cyl, wt, type = "robust")

## -----------------------------------------------------------------------------
# running one-sample proportion test for `vs` at all levels of `am`
mtcars %>%
  group_by(am) %>%
  group_modify(~ contingency_table(.x, vs), .keep = TRUE) %>%
  ungroup()

## -----------------------------------------------------------------------------
knitr::include_graphics("../man/figures/stats_reporting_format.png")

## -----------------------------------------------------------------------------
# needed libraries
library(ggplot2)

# creating a data frame
res <- oneway_anova(iris, Species, Sepal.Length, type = "nonparametric")

ggplot(iris, aes(x = Sepal.Length, y = Species)) +
  geom_boxplot() + # use 'expression' column to display results in the subtitle
  labs(
    x = "Penguin Species",
    y = "Body mass (in grams)",
    title = "Kruskal-Wallis Rank Sum Test",
    subtitle = res$expression[[1]]
  )

