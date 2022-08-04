
#' Refit mixed effects model
#' 
#' Based on a 'GlmmSeq' or 'lmmSeq' class result object, this function attempts
#' to refit an identical model for a specific gene based on the data and fitting
#' parameters stored in the results object and refitting using either
#' [lme4::glmer()] for `GlmmSeq` objects or `lmer()` for `lmmSeq` objects. The
#' fitted model can then be passed on to other packages such as `emmeans` to
#' look at estimated marginal means for the model.
#' 
#' @param object A fitted results object of class `GlmmSeq` or `lmmSeq`
#' @param gene A character value specifying a single gene to extract a fitted
#'   model for
#' @param ... Optional arguments passed to either [lme4::glmer] or [lme4::lmer]
#' @return Fitted model of class `lmerMod` in the case of LMM or `glmerMod` for
#'   a GLMM
#' @export

lmmRefit <- function(object, gene, ...) {
  UseMethod("lmmRefit")
}


#' @export

lmmRefit.lmmSeq <- function(object, gene, ...) {
  data <- object@metadata
  data[, "gene"] <- unlist(object@maindata[gene, ])
  offset <- object@info$offset
  control <- eval(object@info$control)
  fit <- if (object@info$test.stat == "Wald") {
    lme4::lmer(object@formula, data = data,
               control = control, offset = offset, ...)
  } else {
    lmerTest::lmer(object@formula, data = data,
                   control = control, offset = offset, ...)
  }
  fit
}


#' @export

lmmRefit.GlmmSeq <- function(object, gene, ...) {
  data <- object@metadata
  data[, "count"] <- unlist(object@countdata[gene, ])
  disp <- object@info$dispersions[gene]
  offset <- object@info$offset
  control <- eval(object@info$control)
  fit <- lme4::glmer(object@formula, data = data,
                     control = control, offset = offset,
                     family = MASS::negative.binomial(theta = 1/disp),
                     ...)
  fit
}
