#' Glmm for sequencing results of a single gene
#'
#' @param modelFormula the formula
#' @param countdata the sequencing data
#' @param gene the gene name 
#' @param metadata a data frame of sample information
#' @param id Column name in metadata which contains the sample IDs to be used
#' in pairing
#' @param dispersion a numeric for the gene dispersion
#' @param sizeFactor the gene size factor (default=NULL)
#' @param reducedFormula Reduced design formula (default="")
#' @param modelData something something
#' @param control the glmer control (default=glmerControl(optimizer="bobyqa"))
#' @param removeSingles whether to remove unpaired individuals (default=TRUE)
#' @return Returns the fit for the general linear mixed model of a single gene
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl nobars
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom car Anova
#' @importFrom methods slot new
#' @importFrom stats AIC complete.cases logLik reshape terms vcov
#' @keywords hplot
#' @export
#'

glmmGene <- function(modelFormula,
                     countdata,
                     gene,
                     metadata,
                     id,
                     dispersion,
                     sizeFactor=NULL,
                     reducedFormula="",
                     modelData=NULL,
                     control=glmerControl(optimizer="bobyqa"),
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

  
  if (!is.null(sizeFactors)) offset <- log(sf) else offset <- NULL
  cat(paste0('n=', length(ids), ' samples, ', length(unique(ids)), 
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




