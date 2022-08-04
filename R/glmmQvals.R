#' Glmm Sequencing qvalues
#'
#' Add qvalue columns to the glmmSeq dataframe
#' @param object A glmmSeq/lmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}.
#' @param cutoff Prints a table showing the number of probes considered
#' significant by the pvalue cut-off (default=0.05)
#' @param pi0 It is recommended not to input an estimate of pi0. Experienced
#' users can use their own methodology to estimate the proportion of true nulls
#' or set it equal to 1 for the BH procedure (default = NULL).
#' @param verbose Logical whether to print the number of significant probes
#' (default=TRUE)
#' @return Returns a GlmmSeq object with results for gene-wise general linear
#' mixed models with adjusted p-values using the qvalue function
#' @importFrom qvalue qvalue
#' @export
#' @examples
#' data(PEAC_minimal_load)
#' disp <- apply(tpm, 1, function(x) {
#' (var(x, na.rm=TRUE)-mean(x, na.rm = TRUE))/(mean(x, na.rm = TRUE)**2)
#' })
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm[1:5, ],
#'                      metadata = metadata,
#'                      dispersion = disp[1:5],
#'                      verbose=FALSE)
#' MS4A1glmm <- glmmQvals(MS4A1glmm)

glmmQvals <- function(object, cutoff = 0.05, verbose = TRUE) {

  if (!(is(object, "GlmmSeq") | is(object, "lmmSeq"))) {
    stop("object must be a GlmmSeq or lmmSeq object")
  }

  pvals <- object@stats$pvals
  qvals <- apply(pvals, 2, function(x) {
    out <- rep_len(NA, length(x))
    qv <- try(qvalue(x[!is.na(x)])$qvalues, silent = TRUE)
    if (class(qv) != 'try-error') {
      out[!is.na(x)] <- qv
    } else {
      out[!is.na(x)] <- p.adjust(x[!is.na(x)], method='BH')
    }
    out
  })
  rownames(qvals) <- rownames(pvals)
  
  if (verbose) {
    for (i in colnames(qvals)) {
      cat(paste0("\n", i, "\n"))
      cat(paste(rep("-", nchar(i)), collapse = ""))
      print(table(ifelse(qvals[, i] < cutoff,
                         "Significant", "Not Significant")))
    }
  }
  
  object@stats$qvals <- qvals
  return(object)
}
