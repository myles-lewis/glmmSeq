#' Glmm gene-wise
#'
#' Paired plots to show differences between groups and over time
#' @param countdata The count data
#' @param id A character vector of subjects ids for pairing over time.
#' @param time A character vector for the timepoint
#' @param group A character vector for the comparison grouping
#' @param col The group name to be used
#' @param cores The number of cores to use (default=detectCores()/2)
#' @param dispersions A numeric vector of gene dispersions
#' @param removeSingles Logical whether to use only paired samples
#' (default=TRUE)
#' @return Returns dataframe with results for gene-wise glm
#' @importFrom MASS negative.binomial
#' @importFrom lme4 subbars findbars glmer fixef glmerControl
#' @importFrom stats update.formula model.matrix predict setNames
#' @importFrom parallel mclapply
#' @importFrom car Anova
#' @importFrom methods slot
#' @keywords hplot
#' @export
#'

glmmGenes <- function(countdata,
                        id,
                        time,
                        group,
                        col,
                        dispersions,
                        cores=detectCores()/2,
                        removeSingles=TRUE) {


  # Define error messages to catch
  errorA <- paste("Error in (function (fr, X, reTrms, family, nAGQ = 1L,",
                  "verbose = 0L, maxit = 100L, ")
  errorB <- paste("Error in pwrssUpdate(pp, resp, tol = tolPwrss,",
                  "GQmat = GHrule(0L), compDev = compDev, ")
  errorC <- ": \n  pwrssUpdate did not converge in (maxit) iterations\n"
  errorD <- paste(": \n  (maxstephalfit) PIRLS step-halvings failed to reduce",
                  "deviance in pwrssUpdate\n")
  errorE <- ": \n  Downdated VtV is not positive definite\n"

  tryErrors <- c(paste(errorA, errorD),
                 paste(errorA, errorC),
                 paste(errorB, errorC),
                 paste(errorB, errorE),
                 "Error : Response is constant\n", "Unknown error")


  paired <- names(table(id)[table(id) > 1])
  pairedIndex <- as.character(id) %in% paired
  if (removeSingles) {
    countdata <- countdata[, pairedIndex]
    id <- as.character(id)[pairedIndex]
    time <- as.character(time)[pairedIndex]
    group <- as.character(group)[pairedIndex]
  }

  # Check numbers
  if(! all(sapply(list(length(id), length(time)), FUN=identical,
                  length(group)))) stop("Alignment error")
  if(! identical(length(dispersions), nrow(countdata))) {
    stop("Dispersion length must match nrow in countdata")
  }

  # Calculate model
  modelData <- expand.grid(time=unique(time),
                         group=levels(droplevels(factor(group))))
  cat(paste0('n=', length(id), ', ', length(unique(id)), ' individuals\n'))
  designmat <- model.matrix( ~ time * group, modelData)
  start <- Sys.time()

  geneList = lapply(rownames(countdata), function(i) {
    list(y=as.numeric(countdata[i,]), dispersion=dispersions[i])
  })
  cat(paste(length(geneList), "genes"))

  resultList <- mclapply(geneList, function(ylist) {
    data <- data.frame(id=id, time=time, group=group, y=ylist$y)
    fit <- try(
      glmer(y ~ time * group + (1 | id),
            data=data,
            glmerControl(optimizer="bobyqa"),
            family= MASS::negative.binomial(theta=1/ylist$dispersion)),
      silent=F)
    if (class(fit)!='try-error') {
      wald <- car::Anova(fit)
      newY <- predict(fit, newdata=modelData, re.form=NA)
      a <- designmat %*% vcov(fit)
      b <- as.matrix(a %*% t(designmat))
      predvar <- diag(b)
      newSE <- sqrt(predvar)
      newLCI <- exp(newY - newSE * 1.96)
      newUCI <- exp(newY + newSE * 1.96)
      singular <- as.numeric(isSingular(fit))
      convergence <- length(slot(fit, 'optinfo')$conv$lme4$messages)
      return(c(wald[,1], wald[,3], exp(newY),
               newLCI, newUCI, singular, convergence, 0) )
    } else {
      if (fit[1] %in% tryErrors) {errorNumber <- which(tryErrors==fit[1])
      } else {errorNumber <- which(tryErrors=='Unknown error')}
      return( c(rep(NA, 20), try=errorNumber) )
    }
  }, mc.cores=cores)

  # Print timing
  end <- Sys.time()
  print(end - start)

  # Output results
  resultData <- data.frame(t(simplify2array(resultList)))
  rownames(resultData) <- rownames(countdata)
  maxY <- max(as.numeric(droplevels(factor(group))), na.rm=T) * 2
  colnames(resultData) <- c('Wald_time', paste0('Wald_', col),
                          paste0('Wald_', col, ':time'),
                          'P_time', paste0('P_', col),
                          paste0('P_', col, ':time'),
                          paste0('y', 1:maxY), paste0('LCI', 1:maxY),
                          paste0('UCI', 1:maxY),
                          'Singular', 'Convergence.messages', 'Try.error')

  resultData
}


