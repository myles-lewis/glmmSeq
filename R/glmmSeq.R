setClassUnion("character_or_list", c("character", "list"))
setClassUnion("df_or_matrix", c("data.frame", "matrix"))

#' An S4 class to define the glmmSeq output
#'
#' @slot formula The model formula
#' @slot stats the statistics from the glmm fit
#' @slot predict The predicted interception values
#' @slot reducedFormula The reduced formula with removed random effects
#' @slot countdata The input expression data
#' @slot metadata The input metadata
#' @slot modelData the model data for the glmm
#' @slot optInfo Information on whether the model was singular or converged
#' @slot errors Any errors
#' @slot variables The variables used in the formula

setClass("GlmmSeq", slots = list(
  formula = "formula",
  stats = "df_or_matrix",
  predict = "df_or_matrix",
  reducedFormula = "formula",
  countdata = "df_or_matrix",
  metadata = "df_or_matrix",
  modelData = "df_or_matrix",
  optInfo = "matrix",
  errors = "character_or_list",
  variables = "character_or_list"
))


#' Glmm for sequencing results
#'
#' @param modelFormula the model formula. For more information of formula 
#' structure see \code{\link[lme4:glmer]{lme4::glmer()}}  
#' @param countdata the sequencing count data
#' @param metadata a data frame of sample information
#' @param id Column name in metadata which contains the sample IDs to be used
#' in pairing samples
#' @param dispersion a numeric vector of gene dispersion
#' @param sizeFactors size factors (default=NULL). If provided the glmer offset
#' is set to log(sizeFactors). For more information see
#'  \code{\link[lme4:glmer]{lme4::glmer()}}
#' @param reducedFormula Reduced design formula (default="")
#' @param modelData Expanded design matrix
#' @param control the glmer control (default=glmerControl(optimizer="bobyqa")). 
#' For more information see 
#' \code{\link[lme4:glmerControl]{lme4::glmerControl()}}.
#' @param cores number of cores to use. Default=detectCores()/2
#' @param pairedOnly whether to remove unpaired individuals (default=TRUE)
#' @param zeroCount numerical value to offset zeroes for the purpose of log
#' (default=0.125)
#' @param verbose Logical whether to display messaging (default=TRUE)
#' @param ... Other parameters to pass to 
#' \code{\link[lme4:glmer]{lme4::glmer()}}
#' @return Returns a GlmmSeq object with results for gene-wise general linear 
#' mixed models
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars isSingular
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom parallel mclapply detectCores parLapply makeCluster clusterEvalQ 
#' clusterExport stopCluster
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov
#' @export
#' @examples 
#' data(PEAC_minimal_load)
#' 
#' disp <- apply(tpm, 1, function(x){ 
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2) 
#' })
#' 
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      id = 'PATID',
#'                      countdata = tpm['MS4A1', ],
#'                      metadata = metadata,
#'                      dispersion = disp['MS4A1'],
#'                      verbose=FALSE)
#'                      
#' names(attributes(MS4A1glmm))


glmmSeq <- function(modelFormula,
                    countdata,
                    metadata,
                    id,
                    dispersion,
                    sizeFactors=NULL,
                    reducedFormula="",
                    modelData=NULL,
                    control=glmerControl(optimizer="bobyqa"),
                    cores=1,
                    pairedOnly=TRUE,
                    zeroCount=0.125,
                    verbose=TRUE, 
                    ...) {
  
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
  if (zeroCount > 0) countdata[countdata==0] <- zeroCount
  
  # Manipulate formulae
  fullFormula <- update.formula(modelFormula, count ~ ., simplify=FALSE)
  nonRandomFormula <- subbars(modelFormula)
  variables <- rownames(attr(terms(nonRandomFormula), 'factors'))  
  subsetMetadata <- metadata[, variables]  
  ids <- as.character(metadata[, id])
  
  
  # Option to subset to paired samples only
  if (pairedOnly) {
    # Check the distribution for duplicates
    check <- data.frame(table(subsetMetadata))
    check <- check[! check$Freq %in% c(0, 1), ]
    if(nrow(check) > 0){
      mCheck <- as.character(apply(subsetMetadata[, variables], 1, function(x) {
        paste(as.character(x), collapse=" ")
      }))
      cCheck <- as.character(apply(check[, variables], 1, function(x) {
        paste(as.character(x), collapse=" ")
      }))
      countdata <- countdata[, ! mCheck %in% cCheck]
      subsetMetadata <- subsetMetadata[! mCheck %in% cCheck, ]
      ids <- subsetMetadata[, id]
      warning(paste0(paste(check[, id], collapse=", "),
                     " has multiple entries for identical ",
                     paste0(colnames(check)[! colnames(check) %in% 
                                              c(id, "Freq")],
                            collapse=" and "),
                     ". These will be removed."))
    }
    
    paired <- names(table(ids)[table(ids) > 1])
    pairedIndex <- as.character(ids) %in% paired
    countdata <- countdata[, pairedIndex]
    subsetMetadata <- subsetMetadata[pairedIndex, ]
    ids <- ids[pairedIndex]
    sizeFactors <- sizeFactors[pairedIndex]
  }
  
  # Check numbers and alignment
  if(! all(vapply(list(length(ids), nrow(subsetMetadata)), FUN=identical,
                  FUN.VALUE=TRUE, ncol(countdata)))) {
    stop("Alignment error: metadata rownames must match countdata colnames")  
  }
  if(! all(rownames(countdata) %in% names(dispersion), nrow(countdata))) {
    stop("Dispersion length must match nrow in countdata")
  }
  
  if (!is.null(sizeFactors)) offset <- log(sizeFactors) else offset <- NULL
  if (verbose) cat(paste0('\nn=', length(ids), ' samples, ',
                          length(unique(ids)), ' individuals\n'))
  
  
  # setup model prediction
  if (reducedFormula=="") {reducedFormula <- nobars(modelFormula)}
  if (is.null(modelData)) {
    reducedVars <- rownames(attr(terms(reducedFormula), 'factors'))
    varLevels <- lapply(reducedVars, function(x) {
      if (class(metadata[, x])=='factor') {
        return(levels(subsetMetadata[, x]))
      } else {sort(unique(subsetMetadata[, x]))}
    })
    modelData <- expand.grid(varLevels)
    colnames(modelData) <- reducedVars
  }
  designMatrix <- model.matrix(reducedFormula, modelData)
  
  start <- Sys.time()
  fullList <- lapply(rownames(countdata), function(i) {
    list(y=countdata[i,], dispersion=dispersion[i])
  })
  
  # For each gene perform a fit
  if(Sys.info()['sysname'] == "Windows"){
    cl <-makeCluster(cores)
    clusterExport(cl, varlist = c("glmerApply", "fullList", "fullFormula", 
                                  "subsetMetadata", "control", "modelData", 
                                  "offset", "designMatrix", ...), 
                  envir=environment())
    resultList <- parLapply(cl, fullList, function(geneList) {
      glmerApply(geneList, fullFormula = fullFormula, data = subsetMetadata, 
                 control = control, modelData = modelData, offset = offset, 
                 designMatrix = designMatrix, ...)
    })
    stopCluster(cl)
  } else{
    resultList <- mclapply(fullList, function(geneList) {
      glmerApply(geneList, fullFormula=fullFormula, data=subsetMetadata, 
                 control=control, modelData=modelData, offset=offset,
                 designMatrix=designMatrix, ...)
    }, mc.cores=cores)
  }
  
  # Print timing if verbose
  end <- Sys.time()
  if (verbose) print(end - start)
  
  # Output
  names(resultList) <- rownames(countdata)
  noErr <- vapply(resultList, function(x) x$tryErrors=="", FUN.VALUE=TRUE)
  if(length(which(noErr)) == 0) stop("All genes returned an error")
  
  nCheat <- resultList[noErr][[1]]$predict 
  outputPredict <- t(vapply(resultList[noErr], function(x) x$predict,
                            FUN.VALUE=rep(1, length(nCheat))))
  
  outLabels <- apply(modelData, 1, function(x) paste(x, collapse="_"))
  colnames(outputPredict) <- c(paste0("y_", outLabels),
                               paste0("LCI_", outLabels),
                               paste0("UCI_", outLabels))
  
  if (sum(!noErr)!=0) {
    if (verbose) cat(paste0("Errors in ", sum(!noErr), " gene(s):",
                            paste0(names(noErr)[! noErr], collapse=", ")))
    outputErrors <- vapply(resultList[!noErr], function(x) {x$tryErrors},
                           FUN.VALUE = c("test"))
  } else {outputErrors<-c('No errors')}
  
  optInfo <- t(vapply(resultList[noErr], function(x) {
    setNames(x$optinfo, c('Singular', 'Conv'))
  }, FUN.VALUE=c(1,1)))
  
  nCheat <- resultList[noErr][[1]]$stats 
  s <- t(vapply(resultList[noErr], function(x) {x$stats}, 
                FUN.VALUE=rep(1, length(nCheat))))
  
  # Create GlmmSeq object with results
  new("GlmmSeq",
      formula = fullFormula,
      stats = s,
      predict = outputPredict,
      reducedFormula = reducedFormula,
      countdata = countdata,
      metadata = subsetMetadata,
      modelData = modelData,
      optInfo = optInfo,
      errors = outputErrors,
      variables = id
  )
}


#' Fit a glmer model for an individual gene
#'
#' @param geneList List with gene expression and dispersion
#' @param fullFormula the model formula. For more information of formula 
#' structure see \code{\link[lme4:glmer]{lme4::glmer()}}  
#' @param modelData Expanded design matrix
#' @param data The sample data or metadata. 
#' @param designMatrix The design matrix
#' @param control the glmer control (default=glmerControl(optimizer="bobyqa")). 
#' For more information see 
#' \code{\link[lme4:glmerControl]{lme4::glmerControl()}}.
#' @param offset this can be used to specify an a priori known component to be
#'  included in the linear predictor during fitting. For more information see
#'  \code{\link[lme4:glmer]{lme4::glmer()}}. 
#' @param ... Other parameters to pass to 
#' \code{\link[lme4:glmer]{lme4::glmer()}}
#' @return Returns a GlmmSeq object with results for gene-wise general linear 
#' mixed models
#' @importFrom MASS negative.binomial
#' @importFrom lme4 glmer fixef isSingular
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom car Anova
#' @importFrom stats AIC complete.cases logLik reshape terms vcov predict
#' @keywords internal
#' @export
glmerApply <- function(geneList, 
                       fullFormula, 
                       data, 
                       control,
                       modelData, 
                       designMatrix, 
                       offset,
                       ...) {
  data[, 'count'] <- as.numeric(geneList$y)
  fit <- try(suppressMessages(
    lme4::glmer(fullFormula, data=data, control=control, offset=offset,
                family=MASS::negative.binomial(theta=1/geneList$dispersion)), 
    ...),
    silent=TRUE)
  
  if (class(fit)!='try-error') {
    # intercept dropped genes
    if (length(attr(fit@pp$X, "msgRankdrop")) > 0)  {
      return( list(stats=NA, predict=NA, optinfo=NA,
                   tryErrors=attr(fit@pp$X, "msgRankdrop")) )
    }
    stats <- setNames(c(geneList$dispersion, AIC(fit), 
                        as.numeric(logLik(fit))),
                      c('Dispersion', 'AIC', 'logLik'))
    fixedEffects <- lme4::fixef(fit)
    wald <- car::Anova(fit)
    waldtest <- setNames(c(wald[,"Chisq"], wald[,"Pr(>Chisq)"]),
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
    singular <- as.numeric(lme4::isSingular(fit))
    conv <- length(slot(fit, 'optinfo')$conv$lme4$messages)
    rm(fit, data)
    return( list(stats=c(stats, fixedEffects, waldtest),
                 predict=predictdf,
                 optinfo=c(singular, conv),
                 tryErrors="") )
  } else {
    return( list(stats=NA, predict=NA, optinfo=NA, tryErrors=fit[1]) )
  }
}
