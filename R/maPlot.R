#' MA plots
#'
#' @param object A glmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}.
#' @param x1var The name of the first (inner) x parameter
#' @param x2var The name of the second (outer) x parameter
#' @param x1Values Timepoints or categories in `x1var` to be used to calculate fold
#' change. If `NULL` the first two levels in `x1var` are used.
#' @param x2Values Categories in `x2var` to be compared on x and y axis.
#' @param pCutoff The significance cut-off for colour-coding (default=0.01)
#' @param plotCutoff Which probes to include by significance cut-off
#' (default=1 for all markers)
#' @param zeroCountCutoff Which probes to include by minimum counts cut-off
#' (default=50)
#' @param colours Vector of colours to use for significance groups
#' @param labels Row names or indices to label on plot
#' @param fontSize Font size
#' @param labelFontSize Font size for labels
#' @param useAdjusted whether to use adjusted p-values
#' (must have q-values in `object`)
#' @param graphics Either "ggplot" or "plotly"
#' @param verbose Whether to print statistics
#' @return List of three plots. One plot for each `x2Value` and one combined
#' figure
#' @importFrom plotly layout config plot_ly subplot %>%
#' @importFrom ggplot2 ggplot geom_point theme_minimal scale_color_manual labs
#' geom_hline theme element_rect annotate
#' @importFrom ggpubr ggarrange
#' @keywords hplot
#' @export
#' @examples
#' data(PEAC_minimal_load)
#'
#' disp <- apply(tpm, 1, function(x){
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
#' })
#'
#' resultTable <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                        countdata = tpm[1:5, ],
#'                        metadata = metadata,
#'                        dispersion = disp)
#'
#' plots <- maPlot(resultTable,
#'                 x1var='Timepoint',
#'                 x2var='EULAR_6m',
#'                 x2Values=c('Good', 'Non-response'),
#'                 graphics="plotly")
#'
#' plots$combined

maPlot <- function(object,
                   x1var,
                   x2var,
                   x1Values=NULL,
                   x2Values=NULL,
                   pCutoff = 0.01,
                   plotCutoff = 1,
                   zeroCountCutoff = 50,
                   colours=c('grey', 'midnightblue',
                             'mediumvioletred', 'goldenrod'),
                   labels=c(),
                   fontSize=12,
                   labelFontSize=4,
                   useAdjusted=FALSE,
                   graphics="ggplot",
                   verbose=FALSE){

  # Extract the data
  predict <- object@predict
  stats <- object@stats$pvals
  if (useAdjusted) {
    stats <- object@stats$qvals
    colnames(stats) <- paste0("q_", colnames(stats))
  } else {
    stats <- object@stats$pvals
    colnames(stats) <- paste0("P_", colnames(stats))
  }
  adj <- ifelse(useAdjusted, "q_", "P_")
  modelData <- object@modelData
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
    stop("labels must be in rownames(object@predict)")
  }

  # Define the comparisons
  if(is.null(x1Values)){
    x1Values <-  levels(factor(modelData[, x1var]))[1:2]
  }
  if(is.null(x2Values)){
    x2Values <-  levels(factor(modelData[, x2var]))[1:2]
  }
  if(! all(x1Values %in% levels(factor(modelData[, x1var]))) |
     length(x1Values) != 2){
    stop("x1Values must be a vector of two levels in x1var")
  }
  if(! all(x2Values %in% levels(factor(modelData[, x2var]))) |
     length(x2Values) != 2){
    stop("x2Values must be a vector of two levels in x2var")
  }
  xCols <- modelData$y[modelData[, x2var] == x2Values[1] &
                         modelData[, x1var] %in% x1Values]
  yCols <- modelData$y[modelData[, x2var] == x2Values[2] &
                         modelData[, x1var] %in% x1Values]
  plotData$x <- log2(plotData[, xCols[1]]+1) - log2(plotData[, xCols[2]]+1)
  plotData$y <- log2(plotData[, yCols[1]]+1) - log2(plotData[, yCols[2]]+1)
  plotData$maxGroup <- ifelse(abs(plotData$x) > abs(plotData$y),
                              x2Values[1],
                              x2Values[2])
  cols <- gsub("P_", "", colnames(plotData)[grepl("P_", colnames(plotData))])

  # Set up the colour code
  colLevels <- c('Not Significant', paste0(adj, x1var, ' < ', pCutoff),
                 paste0(adj, x2var, " < ", pCutoff, " (up in ",
                        x2Values[2], ")"),
                 paste0(adj, x2var, " < ", pCutoff, " (up in ",
                        x2Values[1], ")"))
  plotData$col <- colLevels[1]
  plotData$col[plotData[, paste0(adj, x1var)] < pCutoff &
                 ! is.na(plotData[, paste0(adj, x1var)])] <- colLevels[2]
  plotData$col[plotData[, paste0(adj, x2var)] < pCutoff &
                 !is.na(plotData[, paste0(adj, x2var)])] <- colLevels[3]
  plotData$col[plotData$col==paste0(adj, x2var, " < ",
                                    pCutoff, " (up in ", x2Values[2], ")") &
                 plotData$maxGroup==x2Values[1]] <- colLevels[4]
  plotData$col[is.na(plotData$col)] <- 'Not Significant'
  plotData$col <- factor(plotData$col, levels=colLevels)

  # Calculate the mean expression
  plotData$meanexp <- apply(
    object@countdata[rownames(object@predict), ], 1, function(x) {
      mean(log2(x+1))
    })
  plotData$zerocount <- apply(
    object@countdata[rownames(object@predict), ], 1, function(x) {
      sum(x < 1)
    })

  # Subset by the plotting cut-off and count cut-off
  plotGenes <- apply(plotData[, grep(adj, colnames(plotData))], 1, function(x){
    any(x < plotCutoff)
  })
  plotData <- plotData[plotGenes, ]
  plotData <- plotData[plotData$zerocount < zeroCountCutoff, ]

  if(verbose){
    cat(paste("Significance\n"))
    print(table(plotData$col))
  }
  labels <- labels[labels %in% rownames(plotData)]

  # Set up annotation for labelled genes
  if (length(labels)!=0) {
    annot1 <- lapply(labels, function(i) {
      row <- plotData[i, ]
      x <- row$meanexp
      y <- row$x
      z <- x^2 + y^2
      list(x=x, y=y,
           text=i, textangle=0, ax=x/z*75, ay=-y/z*75,
           font=list(color="black", size=12),
           arrowcolor="black", arrowwidth=1, arrowhead=0, arrowsize=1.5,
           xanchor="auto", yanchor="auto")
    })
    annot2 <- lapply(labels, function(i) {
      row <- plotData[i, ]
      x <- row$meanexp
      y <- row$y
      z <- x^2 + y^2
      list(x=x, y=y,
           text=i, textangle=0, ax=x/z*75, ay=-y/z*75,
           font=list(color="black", size=labelFontSize),
           arrowcolor="black", arrowwidth=1, arrowhead=0, arrowsize=1.5,
           xanchor="auto", yanchor="auto")
    })
  } else {annot1 <- annot2 <- list()}


  yLab <- paste0("log2(Fold Change)\n", x1var, ": ",
                 x1Values[1], " vs ", x1Values[2])

  # GGplot
  if(graphics == "ggplot"){
    plotData <- plotData[order(plotData$col), ]
    ma1 <- ggplot(data=plotData, aes_string(x="meanexp", y="x", color="col")) +
      geom_point() +
      theme_minimal() +
      scale_color_manual(values=colours, name="") +
      labs(x=bquote(paste("Mean log"[2], "(gene expression + 1)")), y=yLab,
           title=paste0("MA plot (", x2var, " = ", x2Values[1], ")")) +
      geom_hline(yintercept = 0, colour="grey60", linetype="dashed") +
      theme(legend.position=c(1, 0),
            text=element_text(size=fontSize),
            axis.text = element_text(colour = "black", size=fontSize-1),
            legend.background = element_rect(fill=NA, color=NA),
            legend.justification=c(1.1,-0.1)) +
      annotate("text", x=unlist(lapply(annot1, function(x) x$x)),
               y=unlist(lapply(annot1, function(x) x$y)), vjust=1.3,
               size=labelFontSize,
               label= unlist(lapply(annot1, function(x) x$text)))

    ma2 <- ggplot(data=plotData, aes_string(x="meanexp", y="y", color="col")) +
      geom_point() +
      theme_minimal() +
      scale_color_manual(values=colours, breaks=colLevels, name="") +
      labs(x=bquote(paste("Mean log"[2], "(gene expression + 1)")), y=yLab,
           title=paste0("MA plot (", x2var, " = ", x2Values[2], ")")) +
      geom_hline(yintercept = 0, colour="grey60", linetype="dashed") +
      theme(legend.position=c(1, 0),
            text=element_text(size=fontSize),
            axis.text = element_text(colour = "black", size=fontSize-1),
            legend.background = element_rect(fill=NA, color=NA),
            legend.justification=c(1.1, -0.1)) +
      annotate("text", x=unlist(lapply(annot2, function(x) x$x)),
               y=unlist(lapply(annot2, function(x) x$y)), vjust=1.3,
               size=labelFontSize,
               label= unlist(lapply(annot2, function(x) x$text)))

    combined <- ggarrange(ma1, ma2, ncol=1, nrow=2, common.legend = TRUE)

    # Plotly
  } else if(graphics == "plotly"){
    ma1 <- plot_ly(data=plotData, x=~meanexp, y=~x,
                   type='scatter',
                   mode='markers',
                   color=~col,
                   colors=colours,
                   marker=list(size=8, line=list(width=0.5, color='white')),
                   text=rownames(plotData),
                   hoverinfo='text') %>%
      layout(title = paste0("MA plot (", x2var, " = ", x2Values[1], ")"),
             annotations=annot1,
             font=list(size=fontSize),
             xaxis=list(title="Mean log2 gene expr + 1", color='black'),
             yaxis=list(title=yLab, color='black'),
             legend = list(x = 0.65, y = 0.04, font=list(color='black')))

    ma1_dummy <- plot_ly(data=plotData, x=~meanexp, y=~x,
                         type='scatter',
                         mode='markers',
                         color=~col,
                         colors=colours,
                         marker=list(size=8,
                                     line=list(width=0.5, color='white')),
                         text=rownames(plotData),
                         hoverinfo='text', showlegend=FALSE) %>%
      layout(title = paste0("MA plot (", x2var, " = ", x2Values[1], ")"),
             annotations=annot1,
             xaxis=list(title="Mean log2 gene expr + 1", color='black'),
             yaxis=list(title=yLab, color='black'),
             font=list(size=fontSize),
             legend = list(x = 0.65, y = 0.04, font=list(color='black')))

    ma2 <- plot_ly(data=plotData, x=~meanexp, y=~y,
                   type='scatter',
                   mode='markers',
                   color=~col,
                   colors=colours,
                   marker=list(size=8, line=list(width=0.5, color='white')),
                   text=rownames(plotData),
                   hoverinfo='text') %>%
      layout(title = paste0("MA plot (", x2var, " = ", x2Values[2], ")"),
             annotations=annot2,
             font=list(size=fontSize),
             xaxis=list(title="Mean log2 gene expr + 1", color='black'),
             yaxis=list(title=yLab, color='black'),
             legend = list(x = 0.65, y = 0.04, font=list(color='black')))

    combined <- plotly::subplot(ma1_dummy, ma2 %>% layout(title = "MA Plots"),
                                nrows=2,shareY=TRUE,
                                shareX=FALSE,  margin=0.1) %>%
      layout(
        annotations = list(
          list(x = 0, y = 0.45, text = x2Values[1], showarrow = FALSE,
               xref='paper', yref='paper', font=list(size=15)),
          list(x = 0, y = 1.01, text = x2Values[2], showarrow = FALSE,
               xref='paper', yref='paper', font=list(size=15))),
        font=list(size=fontSize),
        legend = list(x = 1, y = 0.5, xanchor="right", yanchor="center")
      )

  } else{stop("graphics must be either ggplot or plotly.")}

  outputList <- list(first = ma1, second = ma2, combined=combined)
  names(outputList)[1:2] <- x2Values

  return(outputList)

}
