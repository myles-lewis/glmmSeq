#' Glmm for sequencing results of a single gene
#'
#' @param modelFormula the formula
#' @param countdata the sequencing data
#' @param gene the gene name
#' @param metadata a data frame of sample information
#' @param id Column name in metadata which contains the sample IDs to be used
#' in pairing
#' @param dispersion a numeric for the gene dispersion
#' @param sizeFactors the gene size factor (default=NULL)
#' @param reducedFormula Reduced design formula (default="")
#' @param modelData something something
#' @param control the glmer control (default=glmerControl(optimizer="bobyqa"))
#' @param zeroCount numerical value to offset zeroes for the purpose of log
#' (default=0.125)
#' @param removeSingles whether to remove unpaired individuals (default=TRUE)
#' @return Returns the fit for the general linear mixed model of a single gene
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov
#' @export
#'

glmmGene <- function(modelFormula,
                     countdata,
                     gene,
                     metadata,
                     id,
                     dispersion,
                     sizeFactors=NULL,
                     reducedFormula="",
                     modelData=NULL,
                     control=glmerControl(optimizer="bobyqa"),
                     zeroCount=0.125,
                     removeSingles=TRUE) {

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
  if(! gene %in% rownames(countdata)){
    stop("gene must be in rownames(countdata)")
  }
  if (! is.numeric(zeroCount)) stop("zeroCount must be numeric")
  if (zeroCount < 0) stop("zeroCount must be >= 0")
  if (zeroCount > 0) countdata[countdata==0] <- zeroCount

  # Manipulate formulae
  fullFormula <- update.formula(modelFormula, count ~ ., simplify=FALSE)
  nf <- subbars(modelFormula)
  variables <- rownames(attr(terms(nf), 'factors'))  # extract variable names
  subsetMetadata <- metadata[, variables]  # restrict metadata to save memory
  ids <- as.character(metadata[, id])
  sf <- sizeFactors

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
    countdata <- countdata[, mCheck != cCheck]
    subsetMetadata <- subsetMetadata[ mCheck != cCheck, ]
    ids <- ids[ids %in% subsetMetadata[, id]]
    warning(paste0(paste(check[, id], collapse=", "),
                   " has multiple entries for identical ",
                   paste0(colnames(check)[! colnames(check) %in% c(id, "Freq")],
                          collapse=" and "),
                   ". These will be removed."))
  }

  # Option to subset to paired samples only
  if (removeSingles) {
    paired <- names(table(ids)[table(ids) > 1])
    pairedIndex <- as.character(ids) %in% paired
    countdata <- countdata[, pairedIndex]
    subsetMetadata <- subsetMetadata[pairedIndex, ]
    ids <- ids[pairedIndex]
    sf <- sizeFactors[pairedIndex]
  }

  # Check numbers and alignment
  if(! all(vapply(list(length(ids), nrow(subsetMetadata)), FUN=identical,
                  FUN.VALUE=TRUE, ncol(countdata)))) {
    stop("Alignment error")
  }


  if (!is.null(sizeFactors)) offset <- log(sf) else offset <- NULL
  cat(paste0('\nn=', length(ids), ' samples, ', length(unique(ids)),
             ' individuals\n'))
  
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

  data <- subsetMetadata
  data[, 'count'] <- as.numeric(countdata[gene,])
  fit <- try(
    glmer(fullFormula, data=data, control=control, offset=offset,
          family=MASS::negative.binomial( theta=1/dispersion)),
    silent=FALSE)

  return(fit)

}




