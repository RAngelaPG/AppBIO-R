## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# library(ggpubr)
# library(agriutilities)
# library(tidyr)
# library(dplyr)
# library(tibble)
# library(asreml)
# 
# head(grassUV) |> print()
# grassUV |>
#   ggplot(
#     aes(x = Time, y = y, group = Plant, color = Plant)
#   ) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~Tmt) +
#   theme_minimal(base_size = 15)

## ----warning=FALSE, message=FALSE, fig.width = 9, echo=FALSE------------------
library(ggpubr)
library(agriutilities)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)

if (requireNamespace("asreml", quietly = TRUE)) {
  library(asreml)
  asreml::asreml.options(trace = 0)
  head(grassUV) |> print()
  grassUV |>
    ggplot(
      aes(x = Time, y = y, group = Plant, color = Plant)
    ) +
    geom_point() +
    geom_line() +
    facet_wrap(~Tmt) +
    theme_minimal(base_size = 15)
}

## ----eval = FALSE-------------------------------------------------------------
# tmp <- grassUV |>
#   group_by(Time, Plant) |>
#   summarise(mean = mean(y, na.rm = TRUE)) |>
#   spread(Time, mean) |>
#   column_to_rownames("Plant")
# 
# a <- covcor_heat(matrix = cor(tmp), legend = "none", size = 4.5) +
#   ggtitle(label = "Pearson Correlation")
# 
# b <- tmp |>
#   cor(use = "pairwise.complete.obs") |>
#   as.data.frame() |>
#   rownames_to_column(var = "Time") |>
#   gather("Time2", "corr", -1) |>
#   type.convert(as.is = FALSE) |>
#   mutate(corr = ifelse(Time < Time2, NA, corr)) |>
#   mutate(Time2 = as.factor(Time2)) |>
#   ggplot(
#     aes(x = Time, y = corr, group = Time2, color = Time2)
#   ) +
#   geom_point() +
#   geom_line() +
#   theme_minimal(base_size = 15) +
#   color_palette(palette = "jco") +
#   labs(color = "Time", y = "Pearson Correlation") +
#   theme(legend.position = "top")
# 
# ggarrange(a, b)

## ----warning=FALSE, message=FALSE, fig.width = 9, echo = FALSE----------------
if (requireNamespace("asreml", quietly = TRUE)) {
  tmp <- grassUV |>
    group_by(Time, Plant) |>
    summarise(mean = mean(y, na.rm = TRUE)) |>
    spread(Time, mean) |>
    column_to_rownames("Plant")

  a <- covcor_heat(matrix = cor(tmp), legend = "none", size = 4.5) +
    ggtitle(label = "Pearson Correlation")

  b <- tmp |>
    cor(use = "pairwise.complete.obs") |>
    as.data.frame() |>
    rownames_to_column(var = "Time") |>
    gather("Time2", "corr", -1) |>
    type.convert(as.is = FALSE) |>
    mutate(corr = ifelse(Time < Time2, NA, corr)) |>
    mutate(Time2 = as.factor(Time2)) |>
    ggplot(
      aes(x = Time, y = corr, group = Time2, color = Time2)
    ) +
    geom_point() +
    geom_line() +
    theme_minimal(base_size = 15) +
    color_palette(palette = "jco") +
    labs(color = "Time", y = "Pearson Correlation") +
    theme(legend.position = "top")
  ggarrange(a, b)
}

## ----eval = FALSE-------------------------------------------------------------
# # Identity variance model.
# model_0 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):idv(Time),
#   data = grassUV
# )
# 
# # Simple correlation model; homogeneous variance form.
# model_1 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):corv(Time),
#   data = grassUV
# )
# 
# # Exponential (or power) model; homogeneous variance form.
# model_2 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):expv(Time),
#   data = grassUV
# )
# 
# # Exponential (or power) model; heterogeneous variance form.
# model_3 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):exph(Time),
#   data = grassUV
# )
# 
# # Antedependence variance model of order 1
# model_4 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):ante(Time),
#   data = grassUV
# )
# 
# # Autoregressive model of order 1; homogeneous variance form.
# model_5 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):ar1v(Time),
#   data = grassUV
# )
# 
# # Autoregressive model of order 1; heterogeneous variance form.
# model_6 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):ar1h(Time),
#   data = grassUV
# )
# 
# # Unstructured variance model.
# model_7 <- asreml(
#   fixed = y ~ Time + Tmt + Tmt:Time,
#   residual = ~ id(Plant):us(Time),
#   data = grassUV
# )

## ----warning=FALSE, message=FALSE, echo=FALSE---------------------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  # Identity variance model.
  model_0 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):idv(Time),
    data = grassUV
  )

  # Simple correlation model; homogeneous variance form.
  model_1 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):corv(Time),
    data = grassUV
  )

  # Exponential (or power) model; homogeneous variance form.
  model_2 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):expv(Time),
    data = grassUV
  )

  # Exponential (or power) model; heterogeneous variance form.
  model_3 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):exph(Time),
    data = grassUV
  )

  # Antedependence variance model of order 1
  model_4 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):ante(Time),
    data = grassUV
  )

  # Autoregressive model of order 1; homogeneous variance form.
  model_5 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):ar1v(Time),
    data = grassUV
  )

  # Autoregressive model of order 1; heterogeneous variance form.
  model_6 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):ar1h(Time),
    data = grassUV
  )

  # Unstructured variance model.
  model_7 <- asreml(
    fixed = y ~ Time + Tmt + Tmt:Time,
    residual = ~ id(Plant):us(Time),
    data = grassUV
  )
}

## ----eval=FALSE---------------------------------------------------------------
# models <- list(
#   "idv" = model_0,
#   "corv" = model_1,
#   "expv" = model_2,
#   "exph" = model_3,
#   "ante" = model_4,
#   "ar1v" = model_5,
#   "ar1h" = model_6,
#   "us" = model_7
# )
# 
# summary_models <- data.frame(
#   model = names(models),
#   aic = unlist(lapply(models, \(x) summary(x)$aic)),
#   bic = unlist(lapply(models, \(x) summary(x)$bic)),
#   loglik = unlist(lapply(models, \(x) summary(x)$loglik)),
#   nedf = unlist(lapply(models, \(x) summary(x)$nedf)),
#   param = unlist(lapply(models, \(x) attr(summary(x)$aic, "param"))),
#   row.names = NULL
# )
# 
# summary_models |> print()
# 
# summary_models |>
#   ggplot(
#     aes(x = reorder(model, -bic), y = bic, group = 1)
#   ) +
#   geom_point(size = 2) +
#   geom_text(aes(x = model, y = bic + 5, label = param), size = 5) +
#   geom_line() +
#   theme_minimal(base_size = 15) +
#   labs(x = NULL, y = "BIC")

## ----warning=FALSE, message=FALSE, fig.width = 9, echo=FALSE------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  models <- list(
    "idv" = model_0,
    "corv" = model_1,
    "expv" = model_2,
    "exph" = model_3,
    "ante" = model_4,
    "ar1v" = model_5,
    "ar1h" = model_6,
    "us" = model_7
  )

  summary_models <- data.frame(
    model = names(models),
    aic = unlist(lapply(models, \(x) summary(x)$aic)),
    bic = unlist(lapply(models, \(x) summary(x)$bic)),
    loglik = unlist(lapply(models, \(x) summary(x)$loglik)),
    nedf = unlist(lapply(models, \(x) summary(x)$nedf)),
    param = unlist(lapply(models, \(x) attr(summary(x)$aic, "param"))),
    row.names = NULL
  )

  summary_models |> print()

  summary_models |>
    ggplot(
      aes(x = reorder(model, -bic), y = bic, group = 1)
    ) +
    geom_point(size = 2) +
    geom_text(aes(x = model, y = bic + 5, label = param), size = 5) +
    geom_line() +
    theme_minimal(base_size = 15) +
    labs(x = NULL, y = "BIC")
}

## ----echo=FALSE---------------------------------------------------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  library(gt)
  models |>
    lapply(\(x) wald(x, denDF = "numeric", ssType = "conditional")$Wald) |>
    lapply(\(x) rownames_to_column(x, "Factor")) |>
    rbindlist(idcol = "Model") |>
    filter(Factor %in% c("Time", "Tmt", "Tmt:Time")) |>
    select(Model, Factor, F.con, Pr) |>
    pivot_wider(names_from = Factor, values_from = c(F.con, Pr)) |>
    mutate_if(is.numeric, round, 3) |>
    gt() |>
    tab_spanner(
      label = "Time",
      columns = c(F.con_Time, Pr_Time)
    ) |>
    tab_spanner(
      label = "Tmt",
      columns = c(F.con_Tmt, Pr_Tmt)
    ) |>
    tab_spanner(
      label = "Tmt:Time",
      columns = c(`F.con_Tmt:Time`, `Pr_Tmt:Time`)
    ) |>
    data_color(
      columns = c(`Pr_Tmt:Time`),
      method = "numeric",
      palette = "ggsci::red_material",
      domain = c(0, 0.08),
      reverse = FALSE
    ) |>
    cols_align(
      align = "center",
      columns = everything()
    ) |>
    cols_label(
      F.con_Tmt = "F.value",
      Pr_Tmt = "P.value",
      `F.con_Tmt:Time` = "F.value",
      `Pr_Tmt:Time` = "P.value",
      F.con_Time = "F.value",
      Pr_Time = "P.value"
    )
}

## ----eval=FALSE---------------------------------------------------------------
# summary(model_4)$varcomp

## ----echo=FALSE---------------------------------------------------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  summary(model_4)$varcomp |> print()
}

## ----eval=FALSE---------------------------------------------------------------
# mat <- extract_rcov(model_4)
# print(mat)

## ----warning=FALSE, message=FALSE, fig.width = 9, echo=FALSE------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  mat <- extract_rcov(model_4)
  print(mat)
}

## ----eval=FALSE---------------------------------------------------------------
# # Plot Correlation  Matrix
# p1 <- covcor_heat(matrix = mat$corr, legend = "none", size = 4.5) +
#   ggtitle(label = "Correlation Matrix (ante)")
# p1
# 
# # Plot Variance-Covariance Matrix
# p2 <- covcor_heat(
#   matrix = mat$vcov,
#   corr = FALSE,
#   legend = "none",
#   size = 4.5,
#   digits = 1
# ) +
#   ggtitle(label = "Covariance Matrix (ante)")
# p2
# ggarrange(p1, p2)

## ----warning=FALSE, message=FALSE, fig.width = 9, echo=FALSE------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  # Plot Correlation  Matrix
  p1 <- covcor_heat(matrix = mat$corr, legend = "none", size = 4.5) +
    ggtitle(label = "Correlation Matrix (ante)")
  p1

  # Plot Variance-Covariance Matrix
  p2 <- covcor_heat(
    matrix = mat$vcov,
    corr = FALSE,
    legend = "none",
    size = 4.5,
    digits = 1
  ) +
    ggtitle(label = "Covariance Matrix (ante)")
  p2
  ggarrange(p1, p2)
}

## ----eval=FALSE---------------------------------------------------------------
# ggarrange(a, p1)

## ----warning=FALSE, message=FALSE, fig.width = 9, echo=FALSE------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  ggarrange(a, p1)
}

## ----eval = FALSE-------------------------------------------------------------
# pvals <- predict(model_4, classify = "Tmt:Time")$pvals
# grassUV |>
#   ggplot(
#     aes(x = Time, y = y, group = Tmt, color = Tmt, shape = Tmt)
#   ) +
#   geom_point(alpha = 0.4, size = 3) +
#   geom_line(data = pvals, mapping = aes(y = predicted.value)) +
#   theme_minimal(base_size = 15) +
#   color_palette(palette = "jco")

## ----warning=FALSE, message=FALSE, fig.width = 9, echo = FALSE----------------
if (requireNamespace("asreml", quietly = TRUE)) {
  pvals <- predict(model_4, classify = "Tmt:Time")$pvals
  grassUV |>
    ggplot(
      aes(x = Time, y = y, group = Tmt, color = Tmt, shape = Tmt)
    ) +
    geom_point(alpha = 0.4, size = 3) +
    geom_line(data = pvals, mapping = aes(y = predicted.value)) +
    theme_minimal(base_size = 15) +
    color_palette(palette = "jco")
}

