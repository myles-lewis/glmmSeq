#' Glmm Sequencing qvalues
#'
#' Add qvalue columns to the glmmSeq dataframe
#' @param glmmResult A glmmSeq object created by 
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
#' 
#' disp <- apply(tpm, 1, function(x){ 
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2) 
#' })
#' 
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      id = 'PATID',
#'                      countdata = tpm[1:10, ],
#'                      metadata = metadata,
#'                      dispersion = disp[1:10],
#'                      verbose=FALSE)
#' 
#' MS4A1glmm <- glmmQvals(MS4A1glmm, pi0=1)                    

glmmQvals <- function(glmmResult, cutoff=0.05, pi0=NULL, verbose=TRUE) {
  
  
  if(class(glmmResult) !="GlmmSeq") stop("glmmResult must be a GlmmSeq object")
  
  resultStats <- data.frame(glmmResult@stats, check.names=FALSE)
  for(cn in colnames(resultStats)[grep('P_', colnames(resultStats))]) {
    q_cn <- gsub('P_', 'q_', cn)
    resultStats[, q_cn] <- NA
    resultStats[!is.na(resultStats[, cn]), q_cn] <-
      qvalue(resultStats[!is.na(resultStats[, cn]), cn], pi0=pi0)$qvalues
    if(verbose){
      cat(paste0("\n", q_cn, "\n"))
      cat(paste(rep("-", nchar(q_cn)), collapse=""))
      print(table(ifelse(resultStats[, q_cn] < cutoff,
                         "Significant", "Not Significant")))
    }
  }
  glmmResult@stats <- as.matrix(resultStats)
  
  return(glmmResult)
}

