#' Plotly Fold Change plot
#'
#' @param glmmResults output from glmmQvals or glmmSeq
#' @param groupCol The grouping column
#' @param colourCutoff The significance cutoff for colour-coding (default=0.01)
#' @param plotCutoff Which probes to include by significance cutoff
#' (default=0.05)
#' @param labels Genes to label on plot
#' @param useAdjust whether to use adjusted pvalues
#' (must have q_ columns in glmmResults)
#' @return Returns dataframe with results for gene-wise glm
#' @importFrom plotly layout config plot_ly
#' @keywords hplot
#' @export
#'

fcPlot <- function(glmmResults,
                   groupCol,
                   colourCutoff=0.01,
                   labels,
                   useAdjusted=F,
                   plotCutoff=0.05){

  predict <- glmmResults@predict
  stats <- glmmResults@stats
  adj <- ifelse(useAdjusted, "q_", "P_")

  plotData <- data.frame(cbind(predict[, c('y1', 'y2', 'y3', 'y4')], stats),
                         check.names = F)

  if(length(grep(adj, colnames(plotData))) < 3) {
    stop(paste("there must be at least 3", adj, "columns"))
  }

  groupVars <- levels(glmmResults@modelData[, groupCol])
  plotData$x <- log2(plotData$y1+1) - log2(plotData$y2+1) # x -> groupVars[1]
  plotData$y <- log2(plotData$y3+1) - log2(plotData$y4+1) # y -> groupVars[2]
  plotData$maxGroup <- ifelse(abs(plotData$x) > abs(plotData$y),
                              groupVars[1],
                              groupVars[2])

  cols <- gsub("P_", "", colnames(plotData)[grepl("P_", colnames(plotData))])
  timeVar <- cols[! grepl(groupCol, cols)]
  interactionVar <- cols[grepl("\\:", cols)]

  plotData$col <- 'Not Significant'
  plotData$col[plotData[, paste0(adj, timeVar)] < colourCutoff &
                 ! is.na(plotData[, paste0(adj, timeVar)])] <-
    paste0(adj, 'time < ', colourCutoff)

  plotData$col[plotData[, paste0(adj, interactionVar)] < colourCutoff &
                 !is.na(plotData[, paste0(adj, interactionVar)])] <-
    paste0(adj, interactionVar, " < ", colourCutoff, " - ", groupVars[2])

  plotData$col[plotData$col==paste0(adj, interactionVar, " < ",
                                    colourCutoff, " - ", groupVars[2]) &
                 plotData$maxGroup==groupVars[1]] <-
    paste0(adj, interactionVar, " < ", colourCutoff, " - ", groupVars[1])
  plotData$col[is.na(plotData$col)] <- 'Not Significant'
  plotData$col <- factor(plotData$col)


  plotGenes <- apply(plotData[, grep(adj, colnames(plotData))], 1, function(x) {
    any(x < plotCutoff)
  })

  print(length(which(plotGenes)))

  plotData <- plotData[plotGenes, ]

  cat(paste("Significance\n"))
  print(table(plotData$col))

  if(any(! labels %in% rownames(plotData))){
    warning(paste(paste(labels[! labels %in% rownames(plotData)], collapse=", "),
            "are not in the glmmResults or do not meet the plotting cutoff so",
            "will not be included in labeling."))
    labels = labels[labels %in% rownames(plotData)]
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




  plot_ly(data=plotData, x=~x, y=~y, type='scatter', mode='markers',
          color=~col, colors=c('grey', 'lightgreen', 'red', 'blue'),
          marker=list(size=8, line=list(width=0.75, color='white')),
          text=rownames(plotData), hoverinfo='text') %>%
    layout(annotations=annot,
           xaxis=list(title=paste0("log2 fold change-", groupVars[1]),
                      color='black'),
           yaxis=list(title=paste0("log2 fold change-", groupVars[2]),
                      color='black'),
           legend = list(x = 0, y = 1, font=list(color='black'))) %>%
    config(edits = list(annotationPosition = TRUE,
                        annotationTail = TRUE,
                        annotationText = TRUE),
           toImageButtonOptions=list(format="svg"))


}
