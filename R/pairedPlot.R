#' Repeated measures plot
#'
#' Plot to show differences between groups over time using base graphics.
#' 
#' @param object A glmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}
#' @param geneName The gene/row name to be plotted
#' @param x1Label The name of the first (inner) x parameter. This must be able
#' to be paired using the ID.
#' @param x2Label The name of the second (outer) x parameter
#' @param IDColumn Column name of sample IDs for pairing
#' @param xTitle Title for the x axis
#' @param yTitle Title for the y axis
#' @param title Plot title. If NULL gene name is used
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shapes (default=19)
#' @param colours The marker colours (default='red') as vector or named vector
#' @param lineColours The line colours (default='grey60') as vector or named 
#' vector
#' @param markerSize Size of markers (default=2)
#' @param fontSize Plot font size
#' @param alpha Line and marker opacity (default=0.7)
#' @param pairedOnly Logical whether to only plot paired samples (default=TRUE)
#' @param addModel Whether to add the fit model with markers (default=TRUE)
#' @param modelSize Size of model points (default=3)
#' @param modelColours Colour of model fit markers (default="black") as vector 
#' or named vector
#' @param modelLineSize Size of model points (default=1) as vector or named 
#' vector
#' @param modelLineColours Colour of model fit lines. If NULL same colours as 
#' modelColour used (default=NULL).
#' @param ... Other parameters to pass to
#' \code{\link[graphics:plot]{graphics::plot()}}
#' @return Returns a paired plot for matched samples
#' @importFrom graphics arrows axis lines mtext plot segments points boxplot
#' @export
#' @examples
#' data(PEAC_minimal_load)
#'
#' disp <- apply(tpm, 1, function(x){
#'   (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
#' })
#'
#' MS4A1glmm <- glmmSeq2(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm['MS4A1', ],
#'                      metadata = metadata,
#'                      dispersion = disp['MS4A1'])
#'
#' pairedPlot(object=MS4A1glmm,
#'            geneName = 'MS4A1',
#'            x1Label = 'Timepoint',
#'            x2Label='EULAR_6m')

pairedPlot <- function(object,
                       geneName = NULL,
                       x1Label = NULL,
                       x2Label = NULL,
                       xTitle = NULL,
                       yTitle = geneName,
                       title=geneName,
                       logTransform=is(object, "GlmmSeq"),
                       shapes=21,
                       colours='grey60',
                       lineColours='grey60',
                       markerSize=0.5,
                       fontSize=NULL,
                       alpha=0.7,
                       pairedOnly=NULL,
                       addModel=TRUE,
                       addPoints=TRUE,
                       modelSize=2,
                       modelColours="royalblue",
                       modelLineSize=1,
                       modelLineColours=modelColours,
                       errorBarLwd=2.5,
                       errorBarLength=0.05,
                       addBox=FALSE,
                       removeDuplicates=FALSE,
                       ...) {
  
  if(!is.numeric(alpha) | alpha > 1 | alpha < 0) {
    stop("alpha must be a number between 0 and 1")
  }
  
  # For outputs from glmmGene
  if(inherits(object, "glmerMod")){
    if(! x1Label %in% colnames(object@frame)) {
      stop("x1Label must be a column name in object@frame")
    }
    if(! x2Label %in% colnames(object@frame)) {
      stop("x2Label must be a column name in object@frame")
    }
    if(! IDColumn %in% colnames(object@frame)) {
      stop("IDColumn must be a column name in object@frame")
    }
    
    x1Values <-  levels(droplevels(factor(object@frame[, x1Label])))
    labs <- levels(droplevels(factor(object@frame[, x2Label])))
    
    # For outputs from glmmSeq
  } else if(inherits(object, "GlmmSeq") |
            inherits(object, "lmmSeq")){
    IDColumn <- object@vars$id
    if (is.null(pairedOnly)) pairedOnly <- object@vars$removeSingles
    if(! x1Label %in% colnames(object@modelData)) {
      stop("x1Label must be a column name in object@modelData")
    }
    if(! x2Label %in% colnames(object@modelData)) {
      stop("x2Label must be a column name in object@modelData")
    }
    maindata <- if (inherits(object, "GlmmSeq")) {
      object@countdata} else object@maindata
    if(! geneName %in% rownames(maindata)){
      stop("geneName must be in rownames(object@countdata)")
    }
    if(ncol(object@modelData) != 2){
      stop(
        paste("These plots only work for interactions between two variable.",
              "Therefore nrow(object@modelData) should be 2."))
    }
    
    # Set up plotting data frame
    df_long <- cbind(
      object@metadata[, c(IDColumn, x1Label, x2Label)],
      geneExp=as.numeric(maindata[geneName, ]))
    pval <- object@stats[geneName, grepl("P_", colnames(object@stats))]
    pval <- formatC(pval, digits=2)
    x1Values <- levels(droplevels(factor(object@modelData[, x1Label])))
    x2Values <- levels(droplevels(factor(object@modelData[, x2Label])))
    labs <- levels(droplevels(factor(object@modelData[, x2Label])))
    df_long$x2 <- as.numeric(factor(df_long[, x2Label]))
    df_long$x1 <- as.numeric(factor(df_long[, x1Label]))
    df_long$x2Labels <- as.character(factor(df_long[, x2Label]))
    df_long$x1Labels <- as.character(factor(df_long[, x1Label]))
    maxX2 <- max(as.numeric(df_long$x2), na.rm=TRUE)
    df_long$x1Factors <- as.numeric(factor(df_long[, x1Label])) +
      (df_long$x2-1) * length(unique(df_long[, x1Label]))
    df_long$id <- df_long[, IDColumn]
    
    # Update colours to named list
    refactorCols <- function(x) {
      setNames(rep(x, length(x2Values))[seq_len(length(x2Values))], x2Values)
    }
    
    # re-factor colours and shapes if necessary
    if(is.null(names(lineColours))) lineColours <- refactorCols(lineColours)
    if(is.null(names(modelColours))) modelColours <- refactorCols(modelColours)
    if(is.null(names(colours))) colours <- refactorCols(colours)
    if(is.null(names(shapes))) shapes <- refactorCols(shapes)
    
    # Set up model fit data
    modelData <- object@modelData
    df_mean <- object@predict[geneName, ]
    outLabels <- apply(modelData, 1, function(x) paste(x, collapse="_"))
    df_mean <- data.frame("x1Factors"=seq_along(modelData[, 1]),
                          "y"= df_mean[ paste0("y_",outLabels)],
                          "lower" =  df_mean[paste0("LCI_", outLabels)],
                          "upper" =  df_mean[paste0("UCI_", outLabels)],
                          "group" = object@modelData[, x2Label])
    
  } else{
    stop("object must be an output from either glmmSeq, lmmSeq or glmmGene")
  }
  
  # Check for duplicated measurements
  checkColumns <-  c(IDColumn, x1Label, x2Label)
  duplicates <- duplicated(df_long[, checkColumns])
  
  if ((removeDuplicates | pairedOnly) & any(duplicates)){
    message("Removing duplicates: multiple entries for ",
            paste(df_long[duplicates, IDColumn], collapse=", "))
    df_long <- df_long[!duplicates, ]
  }
  
  if (length(shapes)==1) shapes <- rep(shapes, maxX2)
  if (length(colours)==1) colours <- rep(colours, maxX2)
  
  if (pairedOnly) {
    # Convert to wide format
    df <- stats::reshape(df_long[, c('id', 'x1', 'geneExp')], timevar='x1',
                         idvar='id', v.names='geneExp', direction='wide')
    # keep only paired samples
    paired <- df$id[complete.cases(df)]
    df_long <- df_long[df_long$id %in% paired, ]
    df_long <- df_long[complete.cases(df_long), ]
  }
  
  if(is.null(title)) title <- geneName
  
  # Generate base plots
  if(logTransform) log <- "y" else log <- ""
  if(is.null(xTitle)) xTitle <- NA
  if(addModel) {
    myYlim <- if (addPoints) {
      range(c(df_mean[, c('lower', 'upper')], df_long$geneExp))
    } else range(df_mean[, c('lower', 'upper')])
  } else myYlim <- NULL
  if (addPoints) {
  plot(as.numeric(df_long$x1Factors), df_long$geneExp,
       ylim = myYlim, type='p', bty='l', las=2,
       xaxt='n', cex.axis=fontSize, cex.lab=fontSize,
       pch=19, 
       col=colours[df_long$x2Labels],
       cex=markerSize, xlab=xTitle, ylab=yTitle,
       log=log,
       ...,
       panel.first={
         for (i in as.character(unique(df_long$id))) {
           lines(df_long$x1Factors[df_long$id==i],
                 df_long$geneExp[df_long$id==i],
                 col=lineColours[unique(df_long$x2Labels[df_long$id == i])] )}
       })
  } else {
    plot(as.numeric(df_long$x1Factors), df_long$geneExp,
         ylim = myYlim, type='n', bty='l', las=2,
         xaxt='n', cex.axis=fontSize, cex.lab=fontSize,
         xlab=xTitle, ylab=yTitle,
         log=log,
         ...)
  }
  if(addModel){
    for(i in as.character(unique(df_mean$group))){
      lines(df_mean$x1Factors[df_mean$group == i],
            df_mean$y[df_mean$group == i],
            lwd=modelSize+1, col=modelLineColours[i])
      arrows(df_mean$x1Factors[df_mean$group == i], 
               df_mean$upper[df_mean$group == i],
               df_mean$x1Factors[df_mean$group == i], 
               df_mean$lower[df_mean$group == i],
               lwd=errorBarLwd, col=modelLineColours[i],
             angle = 90, code = 3, length = errorBarLength)
      points(df_mean$x1Factors[df_mean$group == i], 
             df_mean$y[df_mean$group == i], type = "p",
             col=ifelse(shapes[i] >= 21, "black",
                        modelColours[i]),
             bg=ifelse(shapes[i] < 21, NULL,
                       modelColours[i]),
             pch=shapes[i], cex=modelSize)
    }
  }
  
  if(addBox){
    boxplot(geneExp ~ x1Factors , data = df_long, outcex=0, boxwex=0.25,
            add=TRUE, col=NA, frame = FALSE, axes=FALSE)
  }
  axis(1, seq_along(modelData[, 1]), labels=rep(x1Values, length(x2Values)),
       cex.axis=fontSize)
  axis(1, length(x1Values)*(seq_along(x2Values)-1) +
         length(x1Values)/2+0.5, labels=labs,
       line=1.5, cex.axis=fontSize, tick=FALSE)
  if(title!="") mtext(title, side=3, adj=0, padj=-3, cex=fontSize)
  
  mtext(bquote(
    paste("P"[.(x1Label)]*"=", .(pval[1]),
          ", P"[.(x2Label)]*"=", .(pval[2]),
          ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
          .(pval[3]))), cex=fontSize,
    side=3, adj=0)
  
}
