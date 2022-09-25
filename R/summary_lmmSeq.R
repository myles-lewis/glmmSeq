
#' @export

summary.lmmSeq <- function(object,
                           gene = NULL,
                           digits = max(3L, getOption("digits") - 3L), ...) {
  if (is.null(gene)) {
    statSet <- names(object@stats)
    gp <- lapply(statSet, function(i) {
      out <- object@stats[[i]]
      if (i %in% c("Chisq", "Fval", "Df", "NumDF", "DenDF")) colnames(out) <- paste(i, colnames(out), sep = "_")
      if (i == "pvals") colnames(out) <- paste0("P_", colnames(out))
      if (i == "stdErr") colnames(out) <- paste0("se_", colnames(out))
      out
    })
    do.call(cbind, gp)
  } else {
    out <- lapply(object@stats, function(i) i[gene, ])
    if (is(object, "GlmmSeq")) {
      cat("Generalised linear mixed model\n")
      cat(paste0("Method: ", object@info$method, "\n"))
      if (object@info$method == "lme4") {
        cat("Family: Negative Binomial\n")
      } else {
        cat(paste0("Family: ", object@info$family, "\n"))
      }
    } else {
      cat("Linear mixed model\n")
    }
    print(out$res)
    cat("\nFixed effects:\n")
    cfdf <- data.frame(Estimate = out$coef,
                       `Std. Error` = out$stdErr, check.names = FALSE)
    print(cfdf, digits = digits)
    if (object@info$test.stat == "Wald") {
      cat("\nAnalysis of Deviance Table (Type II Wald chisquare tests)\n")
      testdf <- data.frame(Chisq = out$Chisq,
                           Df = out$Df,
                           `Pr(>Chisq)` = out$pvals, check.names = FALSE)
    } else if (object@info$test.stat == "LRT") {
      cat("\nLikelihood ratio test\n")
      testdf <- data.frame(Chisq = out$Chisq,
                           Df = out$Df,
                           `Pr(>Chisq)` = out$pvals, row.names = " ",
                           check.names = FALSE)
    } else {
      cat("\nType III Analysis of Variance Table with Satterthwaite's method\n")
      testdf <- data.frame(NumDF = out$NumDF,
                           DenDF = out$DenDF,
                           `F value` = out$Fval,
                           `Pr(>F)` = out$pvals, check.names = FALSE)
    }
    print(testdf, digits = digits)
    invisible(out)
  }
}


#' @export

summary.GlmmSeq <- function(object, ...) {
  summary.lmmSeq(object, ...)
}
