#' Plotly or ggplot fold change plots
#'
#' @param object A glmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}.
#' @param x1var The name of the first (inner) x parameter
#' @param x2var The name of the second (outer) x parameter
#' @param x1Values Timepoints or categories in `x1var` used to calculate fold
#' change. If `NULL` the first two levels in `x1var` are used.
#' @param x2Values Categories in `x2var` to be compared on x and y axis.
#' @param pCutoff The significance cut-off for colour-coding
#' (default = 0.01)
#' @param labels Row names or indices to label on plot
#' @param useAdjusted whether to use adjusted p-values (must have q-values in
#'   `object`). Default = FALSE
#' @param plotCutoff Which probes to include on plot by significance cut-off
#' (default = 1, for all markers)
#' @param graphics Graphics system to use: "ggplot" or "plotly"
#' @param fontSize Font size
#' @param labelFontSize Font size for labels
#' @param colours Vector of colours to use for significance groups
#' @param verbose Whether to print statistics
#' @param ... Other parameters to pass to plotly or ggplot
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
#' disp <- apply(tpm, 1, function(x) {
#' (var(x, na.rm = TRUE)-mean(x, na.rm = TRUE))/(mean(x, na.rm = TRUE)**2)
#' })
#'
#' glmmFit <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm[1:5, ],
#'                      metadata = metadata,
#'                      dispersion = disp,
#'                      verbose = FALSE)
#'
#' fcPlot(object = glmmFit,
#'       x1var = "Timepoint",
#'       x2var = "EULAR_6m",
#'       x2Values = c("Good", "Non-response"),
#'       pCutoff = 0.05,
#'       useAdjusted = FALSE,
#'       plotCutoff = 1,
#'       graphics = "plotly")

fcPlot <- function(object,
                   x1var,
                   x2var,
                   x1Values = NULL,
                   x2Values = NULL,
                   pCutoff = 0.01,
                   labels = c(),
                   useAdjusted = FALSE,
                   plotCutoff = 1,
                   graphics = "ggplot",
                   fontSize = 12,
                   labelFontSize = 4,
                   colours = c("grey", "goldenrod1", "red", "blue"),
                   verbose = FALSE, 
                   ...){
  
  # Extract the data
  interactCol <- colnames(object@stats$pvals)
  interactCol <- interactCol[grepl(":", interactCol)]
  if (useAdjusted) {
    stats <- object@stats$qvals
    if (is.null(stats)) stop("Missing q-values", call. = FALSE)
  } else {
    stats <- object@stats$pvals
  }
  adj <- ifelse(useAdjusted, "q_", "P_")
  modelData <- object@modelData

  # Set up the plotting data
  predData <- object@predict[, 1:nrow(modelData)]
  
  if (ncol(stats) < 3) {
    stop("Incorrect p-value table structure")
  }
  if (! all(labels %in% rownames(predData))) {
    stop("Labels not found in rownames(object@predict)")
  }

  # Define the comparisons
  if (is.null(x1Values)) {
    x1Values <- levels(factor(modelData[, x1var]))[1:2]
  }
  if (is.null(x2Values)) {
    x2Values <- levels(factor(modelData[, x2var]))[1:2]
  }
  if (! all(x1Values %in% levels(factor(modelData[, x1var]))) |
     length(x1Values) != 2) {
    stop("x1Values must be a vector of two levels in x1var")
  }
  if (! all(x2Values %in% levels(factor(modelData[, x2var]))) |
     length(x2Values) != 2) {
    stop("x2Values must be a vector of two levels in x2var")
  }
  xCols <- which(modelData[, x2var] == x2Values[1] &
                         modelData[, x1var] %in% x1Values)
  yCols <- which(modelData[, x2var] == x2Values[2] &
                         modelData[, x1var] %in% x1Values)
  
  plotData <- data.frame(
    x = log2(predData[, xCols[2]]+1) - log2(predData[, xCols[1]]+1),
    y = log2(predData[, yCols[2]]+1) - log2(predData[, yCols[1]]+1))
  plotData$maxGroup <- ifelse(abs(plotData$x) > abs(plotData$y),
                              x2Values[1], x2Values[2])

  # Set up the colour code
  colLevels <- c('Not Significant', paste0(adj, x1var, ' < ', pCutoff),
                 paste0(adj, x1var, ":", x2var, " < ", pCutoff,
                        " (biggest FC in ", x2Values[2], ")"),
                 paste0(adj, x1var, ":", x2var, " < ", pCutoff,
                        " (biggest FC in ", x2Values[1], ")"))
  plotData$col <- colLevels[1]
  plotData$col[stats[, x1var] < pCutoff &
                 ! is.na(stats[, x1var])] <- colLevels[2]

  plotData$col[stats[, interactCol] < pCutoff &
                 !is.na(stats[, x2var])] <- colLevels[3]

  plotData$col[plotData$col == colLevels[3] &
                 plotData$maxGroup == x2Values[1]] <- colLevels[4]

  plotData$col[is.na(plotData$col)] <- 'Not Significant'
  plotData$col <- factor(plotData$col, levels = colLevels)

  # Genes passing the cutoff
  plotGenes <- apply(stats, 1, function(x) {
    any(x < plotCutoff)
  })
  plotData <- plotData[plotGenes, ]

  if (verbose) {
    cat("Significance\n")
    print(table(plotData$col))
  }

  # Set up annotations for genes to be labelled
  if (any(! labels %in% rownames(plotData))) {
    warning(
      paste(labels[! labels %in% rownames(plotData)], collapse = ", "),
      "are not in the object or do not meet the plotting cutoff so",
      "will not be included in labeling.")
    labels <- labels[labels %in% rownames(plotData)]
  }
  if (length(labels) != 0) {
    annot <- lapply(labels, function(i) {
      row <- plotData[i, ]
      x <- row$x
      y <- row$y
      z <- sqrt(x^2 + y^2)
      list(x = x, y = y,
           text = i, textangle = 0, ax = x/z*75, ay = -y/z*75,
           font = list(color = "black", size = labelFontSize),
           arrowcolor = "black", arrowwidth = 1, arrowhead = 0, arrowsize = 1.5,
           xanchor = "auto", yanchor = "auto")
    })
  } else {annot <- list()}

  # GGplot
  if (graphics == "ggplot") {
    plotData <- plotData[order(plotData$col), ]
    p <- ggplot(data = plotData, aes_string(x = "x", y = "y", color = "col"), 
                ...) +
      geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+
      geom_point() +
      theme_minimal() +
      scale_color_manual(values = colours, breaks = colLevels, name = "") +
      labs(x = bquote(paste("log"[2], "Fold Change ", .(x1Values[2]), " vs ",
                          .(x1Values[1]), " (", .(x2var), " = ",
                          .(x2Values[1]), ")")),
           y = bquote(paste("log"[2], "Fold Change ", .(x1Values[2]), " vs ",
                          .(x1Values[1]), " (", .(x2var), " = ",
                          .(x2Values[2]), ")")),
           title = "") +
      theme(legend.position = c(0, 1),
            text = element_text(size = fontSize),
            axis.text = element_text(colour = "black", size=fontSize-1),
            legend.background = element_rect(fill = NA, color = NA),
            legend.justification = c(-0.1,1.1),
            plot.margin = unit(c(7, 4, 4, 4), units = "mm")) +
      annotate("text", x = unlist(lapply(annot, function(x) x$x)),
               y = unlist(lapply(annot, function(x) x$y)), vjust = 1.3,
               size = labelFontSize,
               label = unlist(lapply(annot, function(x) x$text)))

    # Plotly
  }else if (graphics == "plotly") {
    plotData <- plotData[order(plotData$col), ]
    p <- plot_ly(data = plotData, x = ~x, y = ~y, type = 'scatter', 
                 mode = 'markers', ...,
                 color = ~col, colors = colours,
                 marker = list(size = 8, 
                               line = list(width = 0.75, color = 'white')),
                 text = rownames(plotData), hoverinfo = 'text') %>%
      layout(annotations = annot,
             xaxis = list(title = paste0("log<sub>2</sub>Fold change ",
                                     x1Values[2], " vs ", x1Values[1],
                                     " (", x2var, " = ", x2Values[1], ")"),
                        color = 'black'),
             yaxis = list(title = paste0("log<sub>2</sub>Fold change ",
                                     x1Values[2], " vs ", x1Values[1],
                                     " (", x2var, " = ", x2Values[2], ")"),
                        color = 'black'),
             font = list(size = fontSize),
             legend = list(x = 0, y = 1, font = list(color = 'black'))) %>%
      config(edits = list(annotationPosition = FALSE,
                          annotationTail = TRUE,
                          annotationText = TRUE),
             toImageButtonOptions = list(format = "svg"))
  } else stop("graphics must be 'ggplot' or 'plotly'")

  return(p)
}
