# Class definitions

setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))
setClassUnion("list_or_matrix", c("list", "matrix"))
setClassUnion("formulaOrNULL", c("formula", "NULL"))

#' An S4 class to define the glmmSeq output
#'
#' @slot info List including the matched call, dispersions, offset, designMatrix
#' @slot formula The model formula
#' @slot stats Statistics from fitted models
#' @slot predict Predicted values
#' @slot reduced Optional reduced formula for LRT
#' @slot countdata The input expression data with count data in rows
#' @slot metadata The input metadata
#' @slot modelData Model data for predictions
#' @slot optInfo Information on whether the model was singular or converged
#' @slot errors Any errors
#' @slot vars List of variables stored from the original call, including the
#'   `id` variable (by default automatically identified from the random effect
#'   term in the model) and `removeSingles` argument

setClass("GlmmSeq", slots = list(
  info = "list",
  formula = "formula",
  stats = "list_or_matrix",
  predict = "df_or_matrix",
  reduced = "formulaOrNULL",
  countdata = "df_or_matrix",
  metadata = "df_or_matrix",
  modelData = "df_or_matrix",
  optInfo = "matrix",
  errors = "character_or_list",
  vars = "list"
))


#' An S4 class to define the lmmSeq output
#'
#' @slot info List including matched call, offset, designMatrix
#' @slot formula The model formula
#' @slot stats Statistics from fitted models
#' @slot predict Predicted values
#' @slot reduced Optional reduced formula for LRT
#' @slot maindata The input expression data with variables in rows
#' @slot metadata The input metadata
#' @slot modelData Model data for predictions
#' @slot optInfo Information on whether the model was singular or converged
#' @slot errors Any errors
#' @slot vars List of variables stored from the original call

setClass("lmmSeq", slots = list(
  info = "list",
  formula = "formula",
  stats = "list_or_matrix",
  predict = "df_or_matrix",
  reduced = "formulaOrNULL",
  maindata = "df_or_matrix",
  metadata = "df_or_matrix",
  modelData = "df_or_matrix",
  optInfo = "matrix",
  errors = "character_or_list",
  vars = "list"
))


ignore_unused_imports <- function() {
  kableExtra::kable_styling()
  emmeans::emmeans()
}