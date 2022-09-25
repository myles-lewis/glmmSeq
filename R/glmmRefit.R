
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
#' @param formula Optional formula to use when refitting model
#' @param control Optional control parameters, see [lme4::lmerControl()] or
#'   [lme4::glmerControl()]
#' @param family Optional GLM family when refitting GLMM using [lme4::glmer()]
#'   or [glmmTMB()]
#' @param ... Optional arguments passed to either [lme4::glmer()] or
#'   [lme4::lmer()]
#' @return Fitted model of class `lmerMod` in the case of LMM, or `glmerMod` or
#'   `glmmTMB` for a GLMM dependent on the original method.
#'   
#' @examples
#' data(PEAC_minimal_load)
#' disp <- apply(tpm, 1, function(x) {
#'   (var(x, na.rm = TRUE)-mean(x, na.rm = TRUE))/(mean(x, na.rm = TRUE)**2)
#' })
#' glmmtest <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm[1:2, ],
#'                      metadata = metadata,
#'                      dispersion = disp,
#'                      verbose = FALSE)
#' 
#' # show summary for single gene
#' summary(glmmtest, "MS4A1")
#' 
#' # refit a single model using lme4::glmer()
#' fit <- glmmRefit(glmmtest, "MS4A1")
#' 
#' # refit model with reduced formula
#' fit2 <- glmmRefit(glmmtest, "MS4A1",
#'                   formula = count ~ Timepoint + EULAR_6m + (1 | PATID))
#' 
#' # LRT
#' anova(fit, fit2)
#' 
#' @export

glmmRefit <- function(object, gene, ...) {
  UseMethod("glmmRefit")
}


#' @rdname glmmRefit
#' @export

glmmRefit.lmmSeq <- function(object, gene,
                             formula = object@formula,
                             ...) {
  data <- object@metadata
  data[, "gene"] <- unlist(object@maindata[gene, ])
  offset <- object@info$offset
  control <- eval(object@info$control)
  fit <- if (object@info$test.stat == "Wald") {
    lme4::lmer(formula, data = data,
               control = control, offset = offset, ...)
  } else {
    lmerTest::lmer(formula, data = data,
                   control = control, offset = offset, ...)
  }
  fit
}


#' @rdname glmmRefit
#' @export

glmmRefit.GlmmSeq <- function(object, gene,
                              formula = object@formula,
                              control = object@info$control,
                              family = NULL,
                              ...) {
  data <- object@metadata
  data[, "count"] <- object@countdata[gene, ]
  offset <- object@info$offset
  if (object@info$method == "lme4") {
    disp <- object@info$dispersion[gene]
    if (is.null(family)) {
      family <- MASS::negative.binomial(theta = 1/disp)
    }
    fit <- lme4::glmer(formula, data = data,
                       control = control, offset = offset,
                       family = family,
                       ...)
  } else {
    if (is.null(family)) {
      family <- eval(object@info$family)
    }
    fit <- glmmTMB(formula, data, family = family,
                   control = control, offset = offset, ...)
  }
  fit
}
