
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
#' @return Fitted model of class `lmerMod` in the case of LMM, or `glmerMod` or `glmmTMB` for
#'   a GLMM dependent on the original method.
#' @export

glmmRefit <- function(object, gene, ...) {
  UseMethod("glmmRefit")
}


#' @export

glmmRefit.lmmSeq <- function(object, gene, ...) {
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

glmmRefit.GlmmSeq <- function(object, gene, ...) {
  data <- object@metadata
  data[, "count"] <- object@countdata[gene, ]
  offset <- object@info$offset
  control <- object@info$control
  if (object@info$method == "lme4") {
    disp <- object@info$dispersion[gene]
    fit <- lme4::glmer(object@formula, data = data,
                       control = control, offset = offset,
                       family = MASS::negative.binomial(theta = 1/disp),
                       ...)
  } else {
    family <- eval(object@info$family)
    fit <- glmmTMB(object@formula, data, family = family,
                   control = control, offset = offset, ...)
  }
  fit
}