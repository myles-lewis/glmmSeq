#' MA plots
#'
#' @param glmmResult output from glmmQvals or glmmSeq
#' @param x1Label The name of the first (inner) x parameter
#' @param x2Label The name of the second (outer) x parameter
#' @param pCutoff The significance cut-off for colour-coding (default=0.01)
#' @param plotCutoff Which probes to include by significance cut-off
#' (default=1 for all markers)
#' @param colours Colour scheme
#' @param labels Genes to label on plot
#' @param transpose Logical whether FC is time or group based (in this case
#' x2Label should be the time parameter)
#' @param useAdjust whether to use adjusted pvalues
#' (must have q_ columns in glmmResult)
#' @param graphics Either "ggplot" or "plotly"
#' @param verbose Whether to print statistics
#' @return List of plot options
#' @importFrom plotly layout config plot_ly subplot %>%
#' @importFrom ggplot2 ggplot geom_point theme_minimal scale_color_manual labs
#' geom_hline theme element_rect annotate
#' @importFrom ggpubr ggarrange
#' @keywords hplot
#' @export
#'

maPlot <- function(glmmResult,
                   x1Label,
                   x2Label,
                   pCutoff = 0.01,
                   plotCutoff = 1,
                   zeroCountCutoff = 50, 
                   colours=c('grey', 'midnightblue', 
                             'mediumvioletred', 'goldenrod'),
                   labels=c(),
                   transpose=FALSE,
                   useAdjusted=FALSE,
                   graphics="ggplot", 
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
  x2Values <- levels(factor(modelData[, x2Label]))
  x1Values <- levels(factor(modelData[, x1Label]))
  plotData$x <- log2(plotData[, xCols[1]]+1) - log2(plotData[, xCols[2]]+1)
  plotData$y <- log2(plotData[, yCols[1]]+1) - log2(plotData[, yCols[2]]+1)
  plotData$maxGroup <- ifelse(abs(plotData$x) > abs(plotData$y),
                              x2Values[1],
                              x2Values[2])
  cols <- gsub("P_", "", colnames(plotData)[grepl("P_", colnames(plotData))])
  
  # Set up the colour code
  plotData$col <- 'Not Significant'
  plotData$col[plotData[, paste0(adj, x1Label)] < pCutoff &
                 ! is.na(plotData[, paste0(adj, x1Label)])] <-
    paste0(adj, x1Label, ' < ', pCutoff)
  plotData$col[plotData[, paste0(adj, x2Label)] < pCutoff &
                 !is.na(plotData[, paste0(adj, x2Label)])] <-
    paste0(adj, x2Label, " < ", pCutoff, " (", x2Values[2], ")")
  plotData$col[plotData$col==paste0(adj, x2Label, " < ",
                                    pCutoff, " (", x2Values[2], ")") &
                 plotData$maxGroup==x2Values[1]] <-
    paste0(adj, x2Label, " < ", pCutoff, " (", x2Values[1], ")")
  plotData$col[is.na(plotData$col)] <- 'Not Significant'
  plotData$col <- factor(plotData$col)
  
  # Calculate the mean expression
  plotData$meanexp <- apply(
    glmmResult@countdata[rownames(glmmResult@predict), ], 1, function(x) {
      mean(log2(x+1))
    })
  plotData$zerocount <- apply(
    glmmResult@countdata[rownames(glmmResult@predict), ], 1, function(x) {
      sum(x < 1)
    })
  
  # Subset by the plotting cutoff and count cutoff
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
           font=list(color="black", size=12),
           arrowcolor="black", arrowwidth=1, arrowhead=0, arrowsize=1.5,
           xanchor="auto", yanchor="auto")
    })
  } else {annot1 <- annot2 <- list()}
  
  
  yLab <- paste0("log2(Fold Change)\n", x1Label, ": ",
                 x1Values[1], " vs ", x1Values[2])
  
  # GGplot
  if(graphics == "ggplot"){
    plotData <- plotData[order(plotData$col), ]
    ma1 <- ggplot(data=plotData, aes_string(x="meanexp", y="x", color="col")) +
      geom_point() +
      theme_minimal() +
      scale_color_manual(values=colours, name="") +
      labs(x=bquote(paste("Mean log"[2], "(gene expression + 1)")), y=yLab,
           title=paste0("MA plot (", x2Label, " = ", x2Values[1], ")")) +
      geom_hline(yintercept = 0, colour="grey60", linetype="dashed") +
      theme(legend.position=c(1, 0), #
            legend.background = element_rect(fill=NA, color=NA),
            legend.justification=c(1.1,-0.1)) +
      annotate("text", x=unlist(lapply(annot1, function(x) x$x)),
               y=unlist(lapply(annot1, function(x) x$y)), vjust=1,
               label= unlist(lapply(annot1, function(x) x$text)))
    
    ma2 <- ggplot(data=plotData, aes_string(x="meanexp", y="y", color="col")) +
      geom_point() +
      theme_minimal() +
      scale_color_manual(values=colours, name="") +
      labs(x=bquote(paste("Mean log"[2], "(gene expression + 1)")), y=yLab,
           title=paste0("MA plot (", x2Label, " = ", x2Values[2], ")")) +
      geom_hline(yintercept = 0, colour="grey60", linetype="dashed") +
      theme(legend.position=c(1, 0),
            legend.background = element_rect(fill=NA, color=NA),
            legend.justification=c(1.1,-0.1)) +
      annotate("text", x=unlist(lapply(annot2, function(x) x$x)),
               y=unlist(lapply(annot2, function(x) x$y)), vjust=1,
               label= unlist(lapply(annot2, function(x) x$text)))
    
    combined <- ggarrange(ma1, ma2, ncol=1, nrow=2, common.legend = T)
    
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
      layout(title = paste0("MA plot (", x2Label, " = ", x2Values[1], ")"), 
             annotations=annot1,
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
      layout(title = paste0("MA plot (", x2Label, " = ", x2Values[1], ")"), 
             annotations=annot1,
             xaxis=list(title="Mean log2 gene expr + 1", color='black'),
             yaxis=list(title=yLab, color='black'),
             legend = list(x = 0.65, y = 0.04, font=list(color='black')))
    
    ma2 <- plot_ly(data=plotData, x=~meanexp, y=~y,
                   type='scatter',
                   mode='markers',
                   color=~col,
                   colors=colours,
                   marker=list(size=8, line=list(width=0.5, color='white')),
                   text=rownames(plotData),
                   hoverinfo='text') %>%
      layout(title = paste0("MA plot (", x2Label, " = ", x2Values[2], ")"), 
             annotations=annot2,
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
        legend = list(x = 1, y = 0.5, xanchor="right", yanchor="center")
      )
    
  } else{stop("graphics must be either ggplot or plotly.")}
  
  outputList <- list(first = ma1, second = ma2, combined=combined)
  names(outputList)[1:2] <- x2Values
  
  return(outputList)
  
}


