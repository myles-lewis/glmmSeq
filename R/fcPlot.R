#' Plotly Fold Change plot
#'
#' @param glmmResult A glmmSeq object created by 
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}.
#' @param x1Label The name of the first (inner) x parameter
#' @param x2Label The name of the second (outer) x parameter
#' @param x1Values to be used to calculate fold change. If NULL the first two 
#' levels in x1Label are used. 
#' @param x2Values Subpopulations in x2Label to be compared on x and y axis.  
#' @param pCutoff The significance cut-off for colour-coding
#' (default=0.01)
#' @param labels Genes to label on plot
#' @param useAdjust whether to use adjusted pvalues
#' (must have q_ columns in glmmResult)
#' @param plotCutoff Which probes to include by significance cut-off
#' (default=1 for all markers)
#' @param graphics Graphics system to use: "ggplot" or "plotly"
#' @param colours Colour vector for significance
#' @param verbose Whether to print statistics
#' @return Returns a plot for fold change between x1Values in one x2Value 
#' subset on x axis and fold change in the other x2Value on the y axis. 
#' @importFrom plotly layout config plot_ly
#' @importFrom ggplot2 ggplot geom_point theme_minimal scale_color_manual labs
#' theme element_rect geom_vline unit
#' @keywords hplot
#' @export
#' @examples 
#' data(PEAC_minimal_load)
#' 
#' disp <- apply(tpm, 1, function(x){ 
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2) 
#' })
#' 
#' glmmFit <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      id = 'PATID',
#'                      countdata = tpm[1:10, ],
#'                      metadata = metadata,
#'                      dispersion = disp,
#'                      verbose=FALSE)
#'                      
#' fcPlot(glmmResult=glmmFit,
#'       x1Label="Timepoint",
#'       x2Label="EULAR_6m",
#'       x2Values=c("Good responder", "Non responder"),
#'       pCutoff=0.05,
#'       useAdjusted = FALSE,
#'       plotCutoff = 1,
#'       graphics="ggplot")

fcPlot <- function(glmmResult,
                   x1Label, 
                   x2Label,
                   x1Values=NULL,
                   x2Values=NULL,
                   pCutoff=0.01,
                   labels=c(),
                   useAdjusted=FALSE,
                   plotCutoff=1, 
                   graphics="ggplot", 
                   colours=c('grey', 'goldenrod1', 'red', 'blue'), 
                   verbose=FALSE){
  
  # Extract the data
  predict <- glmmResult@predict
  stats <- glmmResult@stats
  adj <- ifelse(useAdjusted, "q_", "P_")
  modelData <- glmmResult@modelData
  outLabels <- apply(modelData, 1, function(x) paste(x, collapse="_"))
  modelData$y <- paste0('y_', outLabels)
  
  # Set up the plotting data
  plotData <- data.frame(
    cbind(predict[, paste0('y_', outLabels)], stats),
    check.names = FALSE)
  
  if(length(grep(adj, colnames(plotData))) < 3) {
    stop(paste("there must be at least 3", adj, "columns"))
  }
  if(! all(labels %in% rownames(plotData))){
    stop("labels must be in rownames(glmmResult@predict)")
  }
  
  # Define the comparisons
  if(is.null(x1Values)){
    x1Values <-  levels(factor(modelData[, x1Label]))[1:2]
  } 
  if(is.null(x2Values)){
    x2Values <-  levels(factor(modelData[, x2Label]))[1:2]
  }
  if(! all(x1Values %in% levels(factor(modelData[, x1Label]))) | 
     length(x1Values) != 2){
    stop("x1Values must be a vector of two levels in x1Label")
  } 
  if(! all(x2Values %in% levels(factor(modelData[, x2Label]))) | 
     length(x2Values) != 2){
    stop("x2Values must be a vector of two levels in x2Label")
  } 
  xCols <- modelData$y[modelData[, x2Label] == x2Values[1] & 
                         modelData[, x1Label] %in% x1Values]
  yCols <- modelData$y[modelData[, x2Label] == x2Values[2] & 
                         modelData[, x1Label] %in% x1Values]
  
  plotData$x <- log2(plotData[, xCols[1]]+1) - log2(plotData[, xCols[2]]+1)
  plotData$y <- log2(plotData[, yCols[1]]+1) - log2(plotData[, yCols[2]]+1)
  plotData$maxGroup <- ifelse(abs(plotData$x) > abs(plotData$y),
                              x2Values[1],
                              x2Values[2])
  cols <- gsub("P_", "", colnames(plotData)[grepl("P_", colnames(plotData))])
  cols <- cols[grepl(":", cols)]
  
  # Set up the colour code
  colLevels <- c('Not Significant', paste0(adj, x1Label, ' < ', pCutoff), 
                 paste0(adj, x1Label, ":", x2Label, " < ", pCutoff, 
                        " (biggest FC in ", x2Values[2], ")"), 
                 paste0(adj, x1Label, ":", x2Label, " < ", pCutoff, 
                        " (biggest FC in ", x2Values[1], ")"))
  plotData$col <- colLevels[1]
  plotData$col[plotData[, paste0(adj, cols)] < pCutoff &
                 ! is.na(plotData[, paste0(adj, x1Label)])] <- colLevels[2]
    
  plotData$col[plotData[, paste0(adj,cols)] < pCutoff &
                 !is.na(plotData[, paste0(adj, x2Label)])] <- colLevels[3]
    
  plotData$col[plotData$col==colLevels[3] &
                 plotData$maxGroup==x2Values[1]] <- colLevels[4]
   
  plotData$col[is.na(plotData$col)] <- 'Not Significant'
  plotData$col <- factor(plotData$col, levels=colLevels)
  
  # Genes passing the cutoff
  plotGenes <- apply(plotData[, grep(adj, colnames(plotData))], 1, function(x) {
    any(x < plotCutoff)
  })
  plotData <- plotData[plotGenes, ]
  
  if(verbose){
    cat(paste("Significance\n"))
    print(table(plotData$col))
  }
  
  # Set up annotations for genes to be labelled
  if(any(! labels %in% rownames(plotData))){
    warning(paste(
      paste(labels[! labels %in% rownames(plotData)], collapse=", "),
      "are not in the glmmResult or do not meet the plotting cutoff so",
      "will not be included in labeling."))
    labels <- labels[labels %in% rownames(plotData)]
  }
  if (length(labels)!=0) {
    annot <- lapply(labels, function(i) {
      row <- plotData[i, ]
      x <- row$x
      y <- row$y
      z <- sqrt(x^2 + y^2)
      list(x=x, y=y,
           text=i, textangle=0, ax=x/z*75, ay=-y/z*75,
           font=list(color="black", size=12),
           arrowcolor="black", arrowwidth=1, arrowhead=0, arrowsize=1.5,
           xanchor="auto", yanchor="auto")
    })
  } else {annot <- list()}
  
  # GGplot 
  if(graphics == "ggplot"){
    plotData <- plotData[order(plotData$col), ]
    p <- ggplot(data = plotData, aes_string(x="x", y="y", color="col")) +
      geom_hline(yintercept = 0) + geom_vline(xintercept=0)+
      geom_point() +
      theme_minimal() +
      scale_color_manual(values=colours, breaks = colLevels, name="") +
      labs(x=bquote(paste("log"[2], "Fold Change ", .(x1Values[1]), " vs ", 
                          .(x1Values[2]), " (", .(x2Label), " = ", 
                          .(x2Values[1]), ")")),
           y=bquote(paste("log"[2], "Fold Change ", .(x1Values[1]), " vs ", 
                          .(x1Values[2]), " (", .(x2Label), " = ", 
                          .(x2Values[2]), ")")),
           title="") +
      theme(legend.position=c(0, 1),
            legend.background = element_rect(fill=NA, color=NA),
            legend.justification=c(-0.1,1.1), 
            plot.margin = unit(c(7, 0, 3, 0), units="mm")) +
      annotate("text", x=unlist(lapply(annot, function(x) x$x)),
               y=unlist(lapply(annot, function(x) x$y)), vjust=1,
               label= unlist(lapply(annot, function(x) x$text)))
    
    # Plotly
  }else if(graphics == "plotly"){
    plotData <- plotData[order(plotData$col), ]
    p <- plot_ly(data=plotData, x=~x, y=~y, type='scatter', mode='markers',
                 color=~col, colors=colours,
                 marker=list(size=8, line=list(width=0.75, color='white')),
                 text=rownames(plotData), hoverinfo='text') %>%
      layout(annotations=annot,
             xaxis=list(title=paste0("log<sub>2</sub>Fold change ", 
                                     x1Values[1], " vs ", x1Values[2], 
                                     " (", x2Label, "=", x2Values[1], ")"),
                        color='black'),
             yaxis=list(title=paste0("log<sub>2</sub>Fold change ", 
                                     x1Values[1], " vs ", x1Values[2], 
                                     " (", x2Label, "=", x2Values[2], ")"),
                        color='black'),
             legend = list(x = 0, y = 1, font=list(color='black'))) %>%
      config(edits = list(annotationPosition = TRUE,
                          annotationTail = TRUE,
                          annotationText = TRUE),
             toImageButtonOptions=list(format="svg"))
  } else stop("graphics must be 'ggplot' or 'plotly'")
  
  return(p)
}
