setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the glmmSeq output
#'
#' @slot formula The model formula
#' @slot stats the stats
#' @slot predict The predicted values
#' @slot reducedFormula The reduced formula with removed random effects
#' @slot modelData the model data
#' @slot optinfo the optional info
#' @slot errors errors

setClass("GlmmSeq", slots = list(
  formula = "formula",
  stats = "df_or_matrix",
  predict = "df_or_matrix",
  reducedFormula = "formula",
  countdata = "df_or_matrix",
  modelData = "df_or_matrix",
  optinfo = "matrix",
  errors = "character_or_list"
))


#' Glmm for sequencing results
#'
#' @param modelFormula the formula
#' @param countdata the sequencing data
#' @param metadata a data frame of sample information
#' @param id Column name in metadata which contains the sample IDs to be used
#' in pairing
#' @param dispersions a numeric vector of gene dispersions
#' @param sizeFactors size factors (default=NULL)
#' @param reducedFormula Reduced design formula (default="")
#' @param modelData something something
#' @param control the glmer control (default=glmerControl(optimizer="bobyqa"))
#' @param cores number of cores to use. Default=detectCores()/2
#' @param removeSingles whether to remove unpaired individuals (default=TRUE)
#' @param verbose Logical whether to display messaging (default=TRUE)
#' @return Returns dataframe with results for gene-wise glm
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom parallel mclapply detectCores
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov
#' @keywords hplot
#' @export
#'

glmmSeq <- function(modelFormula,
                    countdata,
                    metadata,
                    id,
                    dispersions,
                    sizeFactors=NULL,
                    reducedFormula="",
                    modelData=NULL,
                    control=glmerControl(optimizer="bobyqa"),
                    cores=detectCores()/2,
                    removeSingles=TRUE,
                    verbose=TRUE) {

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


  # Manipulate formulae
  fullFormula <- update.formula(modelFormula, count ~ ., simplify=FALSE)
  nf <- subbars(modelFormula)
  variables <- rownames(attr(terms(nf), 'factors'))  # extract variable names
  subsetMetadata <- metadata[, variables]  # restrict metadata to save memory
  ids <- as.character(metadata[, id])
  sf <- sizeFactors
  if (removeSingles) {
    paired <- names(table(ids)[table(ids) > 1])
    pairedIndex <- as.character(ids) %in% paired
    countdata <- countdata[, pairedIndex]
    subsetMetadata <- subsetMetadata[pairedIndex, ]
    ids <- ids[pairedIndex]
    sf <- sizeFactors[pairedIndex]
  }

  # Check numbers
  if(! all(vapply(list(length(ids), nrow(subsetMetadata)), FUN=identical,
                  FUN.VALUE=TRUE, ncol(countdata)))) {
    stop("Alignment error")
  }
  if(! identical(length(dispersions), nrow(countdata))) {
    stop("Dispersion length must match nrow in countdata")
  }

  if (!is.null(sizeFactors)) offset <- log(sf) else offset <- NULL
  if (verbose) cat(paste0('n=', length(ids), ' samples, ',
                          length(unique(ids)), ' individuals\n'))


  # setup model prediction
  if (reducedFormula=="") {reducedFormula <- nobars(modelFormula)}
  if (is.null(modelData)) {
    reducedVars <- rownames(attr(terms(reducedFormula), 'factors'))
    varLevels <- lapply(reducedVars, function(x) {
      if (class(metadata[, x])=='factor') return(levels(metadata[, x]))
      sort(unique(metadata[, x]))
    })
    modelData <- expand.grid(varLevels)
    colnames(modelData) <- reducedVars
  }
  designMatrix <- model.matrix(reducedFormula, modelData)

  start <- Sys.time()
  geneList <- lapply(rownames(countdata), function(i) {
    list(y=countdata[i,], dispersion=dispersions[i])
  })

  resultlist <- mclapply(geneList, function(ylist) {
    data <- subsetMetadata
    data[, 'count'] <- as.numeric(ylist$y)
    fit <- try(suppressMessages(
      glmer(fullFormula, data=data, control=control, offset=offset,
            family=MASS::negative.binomial( theta=1/ylist$dispersion))),
      silent=TRUE)

    if (class(fit)!='try-error') {
      # intercept dropped columns
      if (length(attr(fit@pp$X, "msgRankdrop")) > 0)  {
        return( list(stats=NA, predict=NA, optinfo=NA,
                     tryErrors=attr(fit@pp$X, "msgRankdrop")) )
      }
      stats <- setNames(c(ylist$dispersion, AIC(fit), as.numeric(logLik(fit))),
                        c('Dispersion', 'AIC', 'logLik'))
      fixedEffects <- fixef(fit)
      wald <- Anova(fit)
      waldtest <- setNames(c(wald[,1], wald[,3]),
                           c(paste0('Chisq_', rownames(wald)),
                             paste0('P_', rownames(wald))))
      newY <- predict(fit, newdata=modelData, re.form=NA)
      a <- designMatrix %*% vcov(fit)
      b <- as.matrix(a %*% t(designMatrix))
      predVar <- diag(b)
      newSE <- sqrt(predVar)
      newLCI <- exp(newY - newSE * 1.96)
      newUCI <- exp(newY + newSE * 1.96)
      predictdf <- c(exp(newY), newLCI, newUCI)
      singular <- as.numeric(isSingular(fit))
      conv <- length(slot(fit, 'optinfo')$conv$lme4$messages)
      rm(fit, data)
      return( list(stats=c(stats, fixedEffects, waldtest), predict=predictdf,
                   optinfo=c(singular, conv), tryErrors="") )
    } else {
      return( list(stats=NA, predict=NA, optinfo=NA, tryErrors=fit[1]) )
    }
  }, mc.cores=cores)

  # Print timing if verbose
  end <- Sys.time()
  if (verbose) print(end - start)

  # Output
  names(resultlist) <- rownames(countdata)
  noErr <- vapply(resultlist, function(x) x$tryErrors=="", FUN.VALUE=TRUE)
  outputPredict <- t(vapply(resultlist[noErr], function(x) x$predict,
                            FUN.VALUE=rep(1, 12)))
  ny <- nrow(modelData)
  colnames(outputPredict) <- c(paste0('y', 1:ny),
                               paste0('LCI', 1:ny),
                               paste0('UCI', 1:ny))

  if (sum(!noErr)!=0) {
    if (verbose) cat(paste("Errors in", sum(!noErr), "genes\n"))
    outputErrors <- vapply(resultlist[!noErr], function(x) x$try.errors,
                           FUN.VALUE="character")
  } else {outputErrors<-c('No errors')}

  optInfo <- t(vapply(resultlist[noErr], function(x) {
    setNames(x$optinfo, c('Singular', 'Conv'))
  }, FUN.VALUE=c(1,1)))

  s <- t(vapply(resultlist[noErr], function(x) {x$stats}, FUN.VALUE=rep(1, 13)))

  new("GlmmSeq",
      formula = fullFormula,
      stats = s,
      predict = outputPredict,
      reducedFormula = reducedFormula,
      countdata = countdata,
      modelData = modelData,
      optinfo = optInfo,
      errors = outputErrors
  )
}




