#' Plotly Fold Change plot
#'
#' @param glmmResult output from glmmQvals or glmmSeq
#' @param x1Label The name of the first (inner) x parameter
#' @param x2Label The name of the second (outer) x parameter
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
#' @return Returns dataframe with results for gene-wise glm
#' @importFrom plotly layout config plot_ly
#' @importFrom ggplot2 ggplot geom_point theme_minimal scale_color_manual labs
#' theme element_rect geom_vline
#' @keywords hplot
#' @export
#'

fcPlot <- function(glmmResult,
                   x1Label, 
                   x2Label,
                   pCutoff=0.01,
                   labels,
                   useAdjusted=FALSE,
                   plotCutoff=1, 
                   graphics="ggplot", 
                   colours=c('grey', 'lightseagreen', 'red', 'blue'), 
                   verbose=FALSE){
  
  # Extract the data
  predict <- glmmResult@predict
  stats <- glmmResult@stats
  adj <- ifelse(useAdjusted, "q_", "P_")
  modelData <- glmmResult@modelData
  modelData$y <- paste0('y', seq_along(modelData[, 1]))
  
  # Set up the plotting data
  plotData <- data.frame(
    cbind(predict[, paste0('y', seq_along(modelData[, 1]))], stats),
    check.names = FALSE)
  
  if(length(grep(adj, colnames(plotData))) < 3) {
    stop(paste("there must be at least 3", adj, "columns"))
  }
  
  
  # Define the variables
  xCols <- modelData$y[modelData[, x2Label] == 
                         levels(factor(modelData[, x2Label]))[1]]
  yCols <- modelData$y[modelData[, x2Label] == 
                         levels(factor(modelData[, x2Label]))[2]]
  
  groupVars <- levels(factor(modelData[, x2Label]))
  fcVars <- levels(factor(modelData[, x1Label]))
  plotData$x <- log2(plotData[, xCols[1]]+1) - log2(plotData[, xCols[2]]+1)
  plotData$y <- log2(plotData[, yCols[1]]+1) - log2(plotData[, yCols[2]]+1)
  plotData$maxGroup <- ifelse(abs(plotData$x) > abs(plotData$y),
                              groupVars[1],
                              groupVars[2])
  cols <- gsub("P_", "", colnames(plotData)[grepl("P_", colnames(plotData))])
  
  # Set up the colour code
  plotData$col <- 'Not Significant'
  plotData$col[plotData[, paste0(adj, x1Label)] < pCutoff &
                 ! is.na(plotData[, paste0(adj, x1Label)])] <-
    paste0(adj, x1Label, ' < ', pCutoff)
  plotData$col[plotData[, paste0(adj, x2Label)] < pCutoff &
                 !is.na(plotData[, paste0(adj, x2Label)])] <-
    paste0(adj, x2Label, " < ", pCutoff, " (", groupVars[2], ")")
  plotData$col[plotData$col==paste0(adj, x2Label, " < ",
                                    pCutoff, " (", groupVars[2], ")") &
                 plotData$maxGroup==groupVars[1]] <-
    paste0(adj, x2Label, " < ", pCutoff, " (", groupVars[1], ")")
  plotData$col[is.na(plotData$col)] <- 'Not Significant'
  plotData$col <- factor(plotData$col)
  
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
      scale_color_manual(values=colours, name="") +
      labs(x=bquote(paste("log"[2], "Fold Change ", .(fcVars[1]), " vs ", 
                          .(fcVars[2]), " (", .(x2Label), " = ", 
                          .(groupVars[1]), ")")),
           y=bquote(paste("log"[2], "Fold Change ", .(fcVars[1]), " vs ", 
                          .(fcVars[2]), " (", .(x2Label), " = ", 
                          .(groupVars[2]), ")")),
           title="") +
      theme(legend.position=c(0, 1),
            legend.background = element_rect(fill=NA, color=NA),
            legend.justification=c(-0.1,1.1)) +
      annotate("text", x=unlist(lapply(annot, function(x) x$x)),
               y=unlist(lapply(annot, function(x) x$y)), vjust=1,
               label= unlist(lapply(annot, function(x) x$text)))
    
    # Plotly
  }else if(graphics == "plotly"){
    p <- plot_ly(data=plotData, x=~x, y=~y, type='scatter', mode='markers',
                 color=~col, colors=colours,
                 marker=list(size=8, line=list(width=0.75, color='white')),
                 text=rownames(plotData), hoverinfo='text') %>%
      layout(annotations=annot,
             xaxis=list(title=paste0("log<sub>2</sub>Fold change ", 
                                     fcVars[1], " vs ", fcVars[2], 
                                     " (", x2Label, "=", groupVars[1], ")"),
                        color='black'),
             yaxis=list(title=paste0("log<sub>2</sub>Fold change ", 
                                     fcVars[1], " vs ", fcVars[2], 
                                     " (", x2Label, "=", groupVars[2], ")"),
                        color='black'),
             legend = list(x = 0, y = 1, font=list(color='black'))) %>%
      config(edits = list(annotationPosition = TRUE,
                          annotationTail = TRUE,
                          annotationText = TRUE),
             toImageButtonOptions=list(format="svg"))
  } else stop("graphics must be 'ggplot' or 'plotly'")
  
  return(p)
}
