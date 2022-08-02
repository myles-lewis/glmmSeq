

#' @export

lmmRefit <- function(x, ...) {
  UseMethod("lmmRefit")
}



#' @export

lmmRefit.lmmSeq <- function(obj, gene, ...) {
  data <- obj@metadata
  data[, "gene"] <- unlist(obj@maindata[gene, ])
  offset <- obj@info$offset
  control <- eval(obj@info$control)
  fit <- if (obj@info$test.stat == "Wald") {
    lme4::lmer(obj@formula, data = data,
               control = control, offset = offset, ...)
  } else {
    lmerTest::lmer(obj@formula, data = data,
                   control = control, offset = offset, ...)
  }
  fit
}



#' @export

lmmRefit.GlmmSeq <- function(obj, gene, ...) {
  data <- obj@metadata
  data[, "count"] <- unlist(obj@countdata[gene, ])
  disp <- obj@info$dispersions[gene]
  offset <- obj@info$offset
  control <- eval(obj@info$control)
  fit <- lme4::glmer(obj@formula, data = data,
                     control = control, offset = offset,
                     family = MASS::negative.binomial(theta = 1/disp),
                     ...)
  fit
}
