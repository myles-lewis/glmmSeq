#' Mixed model effects plot using ggplot2
#'
#' Plot to show differences between groups and over time using ggplot2.
#' 
#' @param object A glmmSeq/lmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}} or 
#' \code{\link[glmmSeq:lmmSeq]{glmmSeq::lmmSeq()}}
#' @param geneName The gene/row name to be plotted
#' @param x1var The name of the first (inner) x parameter, typically 'time'.
#'   This is anticipated to have different values when matched by ID.
#' @param x2var The name of an optional second (outer) x parameter, which should be a
#'   factor.
#' @param x2shift Amount to shift along x axis for each level of `x2var`. By
#'   default the function will arrange each level of `x2var` side by side.
#' @param xlab Title for the x axis
#' @param ylab Title for the y axis
#' @param plab Optional character vector of labels for p-values. These must
#'   align with column names in `object@stats$pvals`.
#' @param title Plot title. If NULL gene name is used
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shapes (default=19)
#' @param colours The marker colours as vector
#' @param lineColours The line colours (default='grey60') as vector
#' @param markerSize Size of markers (default=1)
#' @param fontSize Plot font size
#' @param alpha Line and marker opacity (default=0.7)
#' @param x2Offset Vertical adjustment to secondary x-axis labels (default=5)
#' @param addPoints Whether to add underlying data points (default=TRUE)
#' @param addModel Whether to add the fit model with markers (default=TRUE)
#' @param modelSize Size of model points (default=4)
#' @param modelColours Colour of model fit markers (default="blue") as vector
#' @param modelLineSize Size of model points (default=1) as vector
#' @param modelLineColours Colour of model fit lines
#' @param addBox Logical whether to add boxplots for mean and IQR
#' @param ... Other parameters to pass to
#' \code{\link[ggplot2:theme]{ggplot2::theme()}}.
#' @return Returns a paired plot for matched samples.
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme geom_boxplot scale_y_continuous
#' scale_x_continuous aes_string annotate margin element_text rel geom_errorbar
#' @importFrom methods is
#' @export
#' @examples
#' data(PEAC_minimal_load)
#'
#' disp <- apply(tpm, 1, function(x){
#'   (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
#' })
#'
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm['MS4A1', , drop = FALSE],
#'                      metadata = metadata,
#'                      dispersion = disp,
#'                      verbose = FALSE)
#'
#' ggmodelPlot(object = MS4A1glmm,
#'            geneName = 'MS4A1',
#'            x1var = 'Timepoint',
#'            x2var = 'EULAR_6m',
#'            colours = c('skyblue', 'goldenrod1', 'mediumvioletred'))

ggmodelPlot <- function(object,
                        geneName = NULL,
                        x1var = NULL,
                        x2var = NULL,
                        x2shift = NULL,
                        xlab = NULL,
                        ylab = geneName,
                        plab = NULL,
                        title = geneName,
                        logTransform = is(object, "GlmmSeq"),
                        shapes = 19,
                        colours = 'grey60',
                        lineColours = 'grey60',
                        markerSize = 1,
                        fontSize = 12,
                        alpha = 0.7,
                        x2Offset = 5,
                        addPoints = TRUE,
                        addModel = TRUE,
                        modelSize = 4,
                        modelColours = "blue",
                        modelLineSize = 1,
                        modelLineColours = modelColours,
                        addBox = FALSE,
                        ...) {
    
  if (!(is(object, "GlmmSeq") | is(object, "lmmSeq"))) {
    stop("object must be an output from glmmSeq or lmmSeq")}
  
  dfs <- formPlot(object, geneName, x1var, x2var, x2shift)
  df_long <- dfs[[1]]
  df_model <- dfs[[2]]
  xdiff <- dfs[[3]]
  x2shift <- dfs[[4]]
  modelData <- object@modelData
  maxX2 <- max(df_long$x2, na.rm = TRUE)
  df_long$x2 <- factor(df_long$x2)
  if (!is.null(x2var)) {
    x2labs <- levels(droplevels(factor(modelData[, x2var])))
  }
  xlim <- range(c(df_long$x, df_model$x), na.rm = TRUE)
  pval <- object@stats$pvals[geneName, , drop = FALSE]
  pval <- formatC(pval, digits=2)
  
  lineColours <- rep_len(lineColours, maxX2)
  modelColours <- rep_len(modelColours, maxX2)
  colours <- rep_len(colours, maxX2)
  shapes <- rep_len(shapes, maxX2)
  
  if (addPoints) {
    p <- ggplot(df_long, 
                aes_string(x="x", y="y", group="id",
                           shape="x2", color="x2")) +
      geom_point(size=markerSize, alpha=alpha, 
                 colour=colours[df_long$x2]) +
      geom_line(alpha=alpha)
  } else {
    p <- ggplot()
  }
  
  p <- p +
    theme_classic() +
    scale_shape_manual(values=shapes) +
    labs(x=xlab, y=ylab, title=title) +
    theme(legend.position = "none", plot.margin=margin(7,4,14+x2Offset,4),
          text=element_text(size=fontSize),
          axis.text = element_text(colour = "black", size=fontSize-1), ...) +
    scale_x_continuous(labels=modelData[, x1var],
                       breaks=if(x2shift < xdiff) {
                         modelData[, x1var]
                       } else {
                         modelData[, x1var] + (as.numeric(modelData[, x2var])-1) * x2shift
                       }
                       ) +
    coord_cartesian(clip = 'off') + 
    scale_color_manual(values=lineColours) 
  
  if(addBox) {
    p <- p + geom_boxplot(mapping=aes_string(x="x", y="y", group="x"),
                          inherit.aes=FALSE,
                          alpha=alpha*0.7, outlier.shape=NA, width=xdiff/6)
  }
  
  if (x2shift >= xdiff) {
    p <- p + geom_text(data=data.frame(label = x2labs,
                                       x = x2shift*(seq_along(x2labs)-1) + xdiff/2,
                                       y = min(c(modelData$LCI, df_long$y),
                                               na.rm=TRUE)),
                       size=rel(4),
                       mapping=aes_string(label="label", x="x", y="y"), hjust = 0.5,
                       nudge_x=0, vjust=x2Offset, inherit.aes=FALSE) +
      theme(axis.title.x = element_text(vjust=-6))
  }
  
  if(addModel){
    p <- p +
      annotate("line", x = df_model$x, y = df_model$y,
               group=df_model$group, size=modelLineSize,
               color=modelLineColours[as.numeric(df_model$group)]) +
      annotate("errorbar", x = df_model$x, y = df_model$y,
               color=modelLineColours[as.numeric(df_model$group)],
               ymin=df_model$lower, ymax=df_model$upper,
               width=xdiff/6, size=modelLineSize) +
      annotate("point", x = df_model$x, y = df_model$y,
               shape=shapes[as.numeric(df_model$group)], size=modelSize,
               color=modelColours[as.numeric(df_model$group)])
  }
  
  if(logTransform) p <- p + scale_y_continuous(trans='log10')
  
  if (is.null(plab)) plab <- colnames(pval)
  ptext <- lapply(1:ncol(pval), function(i) {
    bquote("P" [.(plab[i])] *"="* .(pval[,i]))
  })
  ptext <- bquote(.(paste(unlist(ptext), collapse = '*", "*')))
  p <- p + labs(subtitle = parse(text = ptext))
  
  return(p)
  
}
