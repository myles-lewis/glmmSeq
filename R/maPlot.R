#' MA plots
#'
#' @param glmmResults output from glmmQvals or glmmSeq
#' @param groupCol The grouping column
#' @param colourCutoff The significance cutoff for colour-coding (default=0.01)
#' @param plotCutoff Which probes to include by significance cutoff
#' (default=0.05)
#' @param colours Colour scheme
#' @param labels Genes to label on plot
#' @param useAdjust whether to use adjusted pvalues
#' (must have q_ columns in glmmResults)
#' @return List of plot options
#' @importFrom plotly layout config plot_ly subplot
#' @keywords hplot
#' @export
#'

maPlotly <- function(glmmResults,
                     groupCol,
                     colourCutoff = 0.01,
                     plotCutoff = 0.05,
                     colours=c('grey' , 'orange', 'red', 'blue'),
                     labels=c(),
                     useAdjusted=F){

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

  plotData$meanexp <- apply(
    glmmResults@countdata[rownames(glmmResults@predict), ], 1, function(x) {
      mean(log2(x+1))
      })
  plotData$zerocount <- apply(
    glmmResults@countdata[rownames(glmmResults@predict), ], 1, function(x) {
      sum(x < 1)
      })


  plotGenes <- apply(plotData[, grep(adj, colnames(plotData))], 1, function(x) {
    any(x < plotCutoff)
  })
  plotData <- plotData[plotGenes, ]
  plotData <- plotData[plotData$zerocount < 50, ]


  cat(paste("Significance\n"))
  print(table(plotData$col))
  labels <- labels[labels %in% rownames(plotData)]

  if (length(labels)!=0) {
    annot1 <- lapply(labels, function(i) {
      row <- plotData[i, ]
      x <- row$meanexp
      y <- row$x
      z <- sqrt(x^2 + y^2)
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
      z <- sqrt(x^2 + y^2)
      list(x=x, y=y,
           text=i, textangle=0, ax=x/z*75, ay=-y/z*75,
           font=list(color="black", size=12),
           arrowcolor="black", arrowwidth=1, arrowhead=0, arrowsize=1.5,
           xanchor="auto", yanchor="auto")
    })
  } else {annot1 <- annot2 <- list()}


  ma1 <- plot_ly(data=plotData, x=~meanexp, y=~x,
                 type='scatter',
                 mode='markers',
                 color=~col,
                 colors=colours,
                 marker=list(size=8, line=list(width=0.5, color='white')),
                 text=rownames(plotData),
                 hoverinfo='text',
                 legendgroup = ~col) %>%
    layout(title = groupVars[1], annotations=annot1,
           xaxis=list(title="Mean log2 gene expr + 1", color='black'),
           yaxis=list(title="log2 fold change", color='black'),
           legend = list(x = 0.65, y = 0.04, font=list(color='black')))

  ma2 <- plot_ly(data=plotData, x=~meanexp, y=~y,
                 type='scatter',
                 mode='markers',
                 color=~col,
                 colors=colours,
                 marker=list(size=8, line=list(width=0.5, color='white')),
                 text=rownames(plotData),
                 hoverinfo='text',
                 legendgroup = ~col) %>%
    layout(title = groupVars[2],annotations=annot2,
           xaxis=list(title="Mean log2 gene expr + 1", color='black'),
           yaxis=list(title="log2 fold change", color='black'),
           legend = list(x = 0.65, y = 0.04, font=list(color='black')))

  ma2Dummy <- plot_ly(data=plotData, x=~meanexp, y=~y,
                      type='scatter',
                      mode='markers',
                      color=~col,
                      colors=colours,
                      marker=list(size=8, line=list(width=0.5, color='white')),
                      text=rownames(plotData),
                      hoverinfo='text',
                      legendgroup = ~col,
                      showlegend=F) %>%
    layout(annotations=annot2,
           xaxis=list(title="Mean log2 gene expr + 1", color='black'),
           yaxis=list(title="log2 fold change", color='black'),
           legend = list(x = 0.65, y = 0.04, font=list(color='black')))

  combined <- plotly::subplot(ma1 %>% layout(title=""), ma2Dummy, nrows=2,
                              shareX=F, titleY = T, titleX=T, margin=0.1) %>%
    layout(
      annotations = list(
        list(x = 0.5 , y = 0.45, text = groupVars[1], showarrow = F,
             xref='paper', yref='paper', font=list(size=20)),
        list(x = 0.5 , y = 1, text = groupVars[2], showarrow = F,
             xref='paper', yref='paper', font=list(size=20)))
    )

  outputList <- list(first = ma1, second=ma2, combined=combined)
  names(outputList)[1:2] <- groupVars

  return(outputList)

}


