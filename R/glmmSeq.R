setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the glmmSeq output
#'
#' @slot formula The model formula
#' @slot stats the stats
#' @slot predict The predicted values
#' @slot reducedFormula The reduced formula with removed random effects
#' @slot modelData the modeldata
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
#' @param reducedFormula Reduced design fomula (default="")
#' @param modelData something something
#' @param control the glmer control (default=glmerControl(optimizer="bobyqa"))
#' @param cores number of cores to use. Default=detectCores()/2
#' @param removeSingles whether to remove unpaired individuals (default=TRUE)
#' @param verbose Logical whether to display messaging (default=TRUE)
#' @return Returns dataframe with results for gene-wise glm
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom parallel mclapply
#' @importFrom car Anova
#' @importFrom methods slot
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
                    removeSingles=T,
                    verbose=T) {

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
  fullFormula <- update.formula(modelFormula, count ~ ., simplify=F)  # complete formula
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
  if(! all(sapply(list(length(ids), nrow(subsetMetadata)), FUN=identical,
                  ncol(countdata)))) stop("Alignment error")
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
  noErr <- sapply(resultlist, function(x) x$tryErrors=="")
  outputPredict <- t(sapply(resultlist[noErr], function(x) x$predict))
  ny <- nrow(modelData)
  colnames(outputPredict) <- c(paste0('y', 1:ny),
                               paste0('LCI', 1:ny),
                               paste0('UCI', 1:ny))

  if (sum(!noErr)!=0) {
    if (verbose) cat(paste("Errors in", sum(!noErr), "genes\n"))
    outputErrors <- sapply(resultlist[!noErr], function(x) x$try.errors)
  } else {outputErrors=c('No errors')}

  methods::new("GlmmSeq",
               formula = fullFormula,
               stats = t(sapply(resultlist[noErr], function(x) x$stats)),
               predict = outputPredict,
               reducedFormula = reducedFormula,
               countdata = countdata,
               modelData = modelData,
               optinfo = t(sapply(resultlist[noErr], function(x) {
                 setNames(x$optinfo, c('Singular', 'Conv'))
               })),
               errors = outputErrors
  )
}




