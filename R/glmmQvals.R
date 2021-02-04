#' Glmm Sequencing qvalues
#'
#' Add qvalue columns to the glmmSeq dataframe
#' @param result output from glmmSeq 
#' @param cutoff Prints a table showing the number of probes considered
#' significant by the pvalue cut-off (default=0.05)
#' @param pi0 It is recommended not to input an estimate of pi0. Experienced
#' users can use their own methodology to estimate the proportion of true nulls
#' or set it equal to 1 for the BH procedure (default = NULL).
#' @param verbose Logical whether to print the number of significant probes
#' (default=TRUE)
#' @return Returns dataframe with results for gene-wise glm with qvalue columns
#' @importFrom qvalue qvalue
#' @export

glmmQvals <- function(result, cutoff=0.05, pi0=NULL, verbose=TRUE) {
  
  
  if(class(result) !="GlmmSeq") stop("result must be a GlmmSeq object")
  
  resultStats <- data.frame(result@stats, check.names=FALSE)
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
  result@stats <- as.matrix(resultStats)
  
  return(result)
}

