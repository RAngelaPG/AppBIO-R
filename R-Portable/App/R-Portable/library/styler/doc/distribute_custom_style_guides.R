## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
styler::specify_transformers_drop(
  spaces = list(style_space_around_tilde = "'~'"),
  tokens = list(resolve_semicolon = "';'")
)

