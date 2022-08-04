

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
