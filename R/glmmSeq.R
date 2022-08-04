setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the glmmSeq output
#'
#' @slot info List including the matched call, dispersions, offset, designMatrix
#' @slot formula The model formula
#' @slot stats Statistics from fitted models
#' @slot predict Predicted values
#' @slot reducedFormula The reduced formula with removed random effects
#' @slot countdata The input expression data with count data in rows
#' @slot metadata The input metadata
#' @slot modelData Model data for predictions
#' @slot optInfo Information on whether the model was singular or converged
#' @slot errors Any errors
#' @slot vars List of variables stored from the original call

setClass("GlmmSeq", slots = list(
  info = "list",
  formula = "formula",
  stats = "list",
  predict = "df_or_matrix",
  reducedFormula = "formula",
  countdata = "df_or_matrix",
  metadata = "df_or_matrix",
  modelData = "df_or_matrix",
  optInfo = "matrix",
  errors = "character_or_list",
  vars = "list"
))


#' GLMM with negative binomial distribution for sequencing count data
#'
#' Fits many generalised linear mixed effects models (GLMM) with negative
#' binomial distribution for analysis of overdispersed count data with random
#' effects. Designed for longitudinal analysis of RNA-Sequencing count data.
#' Wald type 2 Chi-squared test is used to calculate p-values.
#'
#' @param modelFormula the model formula. This must be of the form `"~ ..."`
#'   where the structure is assumed to be `"counts ~ ..."`. The formula must
#'   include a random effects term. For more information on formula structure
#'   for random effects see \code{\link[lme4:glmer]{lme4::glmer()}}
#' @param countdata the sequencing count data matrix with genes in rows and
#'   samples in columns
#' @param metadata a dataframe of sample information with variables in columns
#'   and samples in rows
#' @param id Optional. Used to specify the column in metadata which contains the
#'   sample IDs to be used in repeated samples for random effects. If not
#'   specified, the function defaults to using the variable after the "|" in the
#'   random effects term in the formula.
#' @param dispersion a numeric vector of gene dispersion
#' @param sizeFactors size factors (default = NULL). If provided the `glmer` 
#' offset is set to log(sizeFactors). For more information see``
#'  \code{\link[lme4:glmer]{lme4::glmer()}}
#' @param reducedFormula Reduced design formula (default = "")
#' @param modelData Expanded design matrix
#' @param designMatrix custom design matrix
#' @param control the `glmer` optimizer control (default =
#'   `glmerControl(optimizer = "bobyqa")`). See
#'   \code{\link[lme4:glmerControl]{lme4::glmerControl()}}.
#' @param cores number of cores to use. Default = 1.
#' @param removeSingles whether to remove individuals without repeated measures
#' (default = FALSE)
#' @param zeroCount numerical value to offset zeroes for the purpose of log
#' (default = 0.125)
#' @param verbose Logical whether to display messaging (default = TRUE)
#' @param returnList Logical whether to return results as a list or `glmmSeq` 
#' object (default = FALSE). Useful for debugging.
#' @param progress Logical whether to display a progress bar
#' @param ... Other parameters to pass to
#' \code{\link[lme4:glmer]{lme4::glmer()}}
#' @return Returns an S4 class `GlmmSeq` object with results for gene-wise
#'   general linear mixed models. A list of results is returned if `returnList`
#'   is `TRUE` which is useful for debugging.
#' @details
#' This function is a wrapper for [lme4::glmer()]. Wald type 2 Chi-squared test
#' is calculated as per [car::Anova()] optimised for speed. Parallelisation is
#' provided using [parallel::mclapply] on Unix/Mac or [parallel::parLapply] on
#' PC.
#' 
#' @examples
#' data(PEAC_minimal_load)
#' disp <- apply(tpm, 1, function(x) {
#' (var(x, na.rm = TRUE)-mean(x, na.rm = TRUE))/(mean(x, na.rm = TRUE)**2)
#' })
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm["MS4A1", ],
#'                      metadata = metadata,
#'                      dispersion = disp["MS4A1"],
#'                      verbose = FALSE)
#' names(attributes(MS4A1glmm))
#' 
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars isSingular
#' @importFrom parallel mclapply detectCores parLapply makeCluster clusterEvalQ
#' clusterExport stopCluster
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov pchisq
#'   update.formula model.matrix predict setNames coef
#' @export


glmmSeq <- function(modelFormula,
                    countdata,
                    metadata,
                    id = NULL,
                    dispersion,
                    sizeFactors = NULL,
                    reducedFormula = "",
                    modelData = NULL,
                    designMatrix = NULL,
                    control = glmerControl(optimizer = "bobyqa"),
                    cores = 1,
                    removeSingles = FALSE,
                    zeroCount = 0.125,
                    verbose = TRUE,
                    returnList = FALSE, 
                    progress = FALSE,
                    ...) {
  glmmcall <- match.call(expand.dots = TRUE)
  # Catch errors
  if (length(findbars(modelFormula)) == 0) {
    stop("No random effects terms specified in formula")
  }
  if (ncol(countdata) != nrow(metadata)) {
    stop("countdata columns different size to metadata rows")
  }
  if (!is.null(sizeFactors) & ncol(countdata) != length(sizeFactors)) {
    stop("Different sizeFactors length")
  }
  if (! is.numeric(zeroCount)) stop("zeroCount must be numeric")
  if (zeroCount < 0) stop("zeroCount must be >= 0")
  if (zeroCount > 0) countdata[countdata == 0] <- zeroCount
  
  # Manipulate formulae
  fullFormula <- update.formula(modelFormula, count ~ ., simplify = FALSE)
  nonRandomFormula <- subbars(modelFormula)
  variables <- rownames(attr(terms(nonRandomFormula), "factors"))
  subsetMetadata <- metadata[, variables]
  if (is.null(id)) {
    fb <- findbars(modelFormula)
    id <- sub(".*[|]", "", fb)
    id <- gsub(" ", "", id)
  }
  ids <- as.character(metadata[, id])

  # Option to subset to remove unpaired samples
  if (removeSingles) {
    nonSingles <- names(table(ids))[table(ids) > 1]
    nonSingleIDs <- ids %in% nonSingles
    countdata <- countdata[, nonSingleIDs]
    sizeFactors <- sizeFactors[nonSingleIDs]
    subsetMetadata <- subsetMetadata[nonSingleIDs, ]
    ids <- ids[nonSingleIDs]
  }
  
  if (! all(rownames(countdata) %in% names(dispersion))) {
    stop("Some dispersion values are missing")
  }
  
  if (!is.null(sizeFactors)) offset <- log(sizeFactors) else offset <- NULL
  if (verbose) cat(paste0("\nn = ", length(ids), " samples, ",
                          length(unique(ids)), " individuals\n"))
  
  
  # setup model prediction
  if (reducedFormula == "") reducedFormula <- nobars(modelFormula)
  if (is.null(modelData)) {
    reducedVars <- rownames(attr(terms(reducedFormula), "factors"))
    varLevels <- lapply(reducedVars, function(x) {
      if (is.factor(metadata[, x])) {
        return(levels(subsetMetadata[, x]))
      } else {sort(unique(subsetMetadata[, x]))}
    })
    modelData <- expand.grid(varLevels)
    colnames(modelData) <- reducedVars
  } 

  if (is.null(designMatrix)){
    designMatrix <- model.matrix(reducedFormula, modelData)
  } 
  
  # Adapted from car:::Anova.II.mer
  reduced2 <- nobars(fullFormula)
  fac <- attr(terms(reduced2), "factors")
  data2 <- metadata
  data2[,'count'] <- rep(0, nrow(data2)) 
  dm2 <- model.matrix(reduced2, data2)
  assign <- attr(dm2, "assign")
  term.labels <- attr(terms(reduced2), "term.labels")
  p <- length(assign)
  I.p <- diag(p)
  n.terms <- length(term.labels)
  hyp.matrix.1 <- hyp.matrix.2 <- list()
  for (i in seq_len(n.terms)) {
    which.term <- i
    subs.term <- which(assign == which.term)
    relatives <- car_relatives(term.labels[i], term.labels, fac)
    subs.relatives <- NULL
    for (relative in relatives) subs.relatives <- c(subs.relatives, 
                                                    which(assign == relative))
    hyp.matrix.1[[i]] <- I.p[subs.relatives, , drop = FALSE]
    hyp.matrix.2[[i]] <- I.p[c(subs.relatives, subs.term), , drop = FALSE]
  }
  names(hyp.matrix.1) <- term.labels
  
  start <- Sys.time()
  fullList <- lapply(rownames(countdata), function(i) {
    list(y = as.numeric(countdata[i, ]), dispersion = dispersion[i])
  })
  
  # For each gene perform a fit
  if (Sys.info()["sysname"] == "Windows" & cores > 1) {
    cl <- makeCluster(cores)
    clusterExport(cl, varlist = c("glmerCore", "fullList", "fullFormula",
                                  "subsetMetadata", "control", "modelData",
                                  "offset", "designMatrix",
                                  "hyp.matrix.1", "hyp.matrix.2", ...),
                  envir = environment())
    if (progress) {
      resultList <- pblapply(fullList, function(geneList) {
        glmerCore(geneList, fullFormula = fullFormula, data = subsetMetadata,
                  control = control, modelData = modelData, offset = offset,
                  designMatrix = designMatrix,
                  hyp.matrix.1 = hyp.matrix.1,
                  hyp.matrix.2 = hyp.matrix.2, ...)
      }, cl = cl)
    } else {
      resultList <- parLapply(cl = cl, fullList, function(geneList) {
        glmerCore(geneList, fullFormula = fullFormula, data = subsetMetadata,
                  control = control, modelData = modelData, offset = offset,
                  designMatrix = designMatrix,
                  hyp.matrix.1 = hyp.matrix.1,
                  hyp.matrix.2 = hyp.matrix.2, ...)
      })
    }
    stopCluster(cl)
  } else{
    if (progress) {
      resultList <- pbmclapply(fullList, function(geneList) {
        glmerCore(geneList, fullFormula = fullFormula, data = subsetMetadata,
                  control = control, modelData = modelData, offset = offset,
                  designMatrix = designMatrix,
                  hyp.matrix.1 = hyp.matrix.1,
                  hyp.matrix.2 = hyp.matrix.2, ...)
      }, mc.cores = cores)
      if ("value" %in% names(resultList)) resultList <- resultList$value
    } else {
      resultList <- mclapply(fullList, function(geneList) {
        glmerCore(geneList, fullFormula = fullFormula, data = subsetMetadata,
                  control = control, modelData = modelData, offset = offset,
                  designMatrix = designMatrix,
                  hyp.matrix.1 = hyp.matrix.1,
                  hyp.matrix.2 = hyp.matrix.2, ...)
      }, mc.cores = cores)
    }
  }
  if(returnList) return(resultList)
  
  # Print timing if verbose
  end <- Sys.time()
  if (verbose) print(end - start)
  
  # Output
  names(resultList) <- rownames(countdata)
  noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
  if (length(which(noErr)) == 0) { 
    stop("All genes returned an error. Check sufficient data in each group")
  }
  
  predList <- lapply(resultList[noErr], "[[", "predict")
  outputPredict <- do.call(rbind, predList)
  
  outLabels <- apply(modelData, 1, function(x) paste(x, collapse = "_"))
  colnames(outputPredict) <- c(paste0("y_", outLabels),
                               paste0("LCI_", outLabels),
                               paste0("UCI_", outLabels))
  
  if (sum(!noErr) != 0) {
    if (verbose) cat(paste0("Errors in ", sum(!noErr), " gene(s):",
                            paste0(names(noErr)[! noErr], collapse = ", ")))
    outputErrors <- vapply(resultList[!noErr], function(x) {x$tryErrors},
                           FUN.VALUE = c("test"))
  } else {outputErrors <- c("No errors")}
  
  optInfo <- t(vapply(resultList[noErr], function(x) {
    setNames(x$optinfo, c("Singular", "Conv"))
  }, FUN.VALUE = c(1, 1)))
  
  s <- organiseStats(resultList[noErr], "Wald")
  
  # Create GlmmSeq object with results
  new("GlmmSeq",
      info = list(call = glmmcall,
                  offset = offset,
                  designMatrix = designMatrix,
                  control = substitute(control),
                  test.stat = "Wald",
                  dispersion = dispersion),
      formula = fullFormula,
      stats = s,
      predict = outputPredict,
      reducedFormula = reducedFormula,
      countdata = countdata,
      metadata = subsetMetadata,
      modelData = modelData,
      optInfo = optInfo,
      errors = outputErrors,
      vars = list(id = id,
                  removeSingles = removeSingles)
  )
}


glmerCore <- function(geneList,
                      fullFormula,
                      data,
                      control,
                      modelData,
                      designMatrix,
                      offset,
                      hyp.matrix.1,
                      hyp.matrix.2,
                      ...) {
  data[, "count"] <- geneList$y
  fit <- try(suppressMessages(suppressWarnings(
    lme4::glmer(fullFormula, data = data, control = control, offset = offset,
                family = MASS::negative.binomial(theta = 
                                                   1/geneList$dispersion),
                ...))),
    silent = TRUE)
  
  if (!inherits(fit, "try-error")) {
    # intercept dropped genes
    if (length(attr(fit@pp$X, "msgRankdrop")) > 0)  {
      return( list(stats = NA, predict = NA, optinfo = NA,
                   tryErrors = attr(fit@pp$X, "msgRankdrop")) )
    }
    stats <- setNames(c(geneList$dispersion, AIC(fit),
                        as.numeric(logLik(fit))),
                      c("Dispersion", "AIC", "logLik"))
    fixedEffects <- lme4::fixef(fit)
    stdErr <- suppressWarnings(coef(summary(fit))[, 2])
    vcov. <- suppressWarnings(vcov(fit, complete = FALSE))
    vcov. <- as.matrix(vcov.)
    waldtest <- lmer_wald(fixedEffects, hyp.matrix.1, hyp.matrix.2, vcov.)
    
    newY <- predict(fit, newdata = modelData, re.form = NA)
    a <- designMatrix %*% vcov.
    b <- as.matrix(a %*% t(designMatrix))
    predVar <- diag(b)
    newSE <- sqrt(predVar)
    newLCI <- exp(newY - newSE * 1.96)
    newUCI <- exp(newY + newSE * 1.96)
    predictdf <- c(exp(newY), newLCI, newUCI)
    singular <- as.numeric(lme4::isSingular(fit))
    conv <- length(slot(fit, "optinfo")$conv$lme4$messages)
    rm(fit, data)
    return(list(stats = stats,
                coef = fixedEffects,
                stdErr = stdErr,
                chisq = waldtest$chisq,
                df = waldtest$df,
                predict = predictdf,
                optinfo = c(singular, conv),
                tryErrors = "") )
  } else {
    return(list(stats = NA, coef = NA, stdErr = NA, chisq = NA, df = NA,
                predict = NA, optinfo = NA, tryErrors = fit[1]))
  }
}


#' @export

summary.GlmmSeq <- function(object, ...) {
  summary.lmmSeq(object, ...)
}

