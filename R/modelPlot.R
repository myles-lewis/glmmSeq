#' Model plot
#'
#' Model plots to show the overall differences between groups and over time
#' @param glmmResult A glmmSeq object created by 
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}.. 
#' @param geneName Gene name to plot
#' @param x1Label The name of the first (inner) x parameter
#' @param x2Label The name of the second (outer) x parameter
#' @param xTitle Title for the x axis
#' @param yTitle Title for the y axis
#' @param title Plot title. If NULL gene name is used. 
#' @param markerSize Size of markers (default=5)
#' @param fontSize Plot font size (default=10)
#' @param lineWidth Plot line size (default=1)
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shape (default=21)
#' @param x2Offset Vertical adjustment to secondary x-axis (default=3)
#' @param colours The marker colour, default=c('red', 'blue')
#' @param overlap Logical whether x2Label fits should be plotted overlapping one
#' another.
#' @param addErrorbars Logical whether to add error bars. 
#' @param graphics Which graphic system to use (options = "base" or "ggplot")
#' @param ... Other parameters to pass to 
#' \code{\link[graphics:plot]{graphics::plot()}} or 
#' \code{\link[ggplot2:theme]{ggplot2::theme()}}. 
#' @return Returns a plot with the glmm fit for a given gene
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme scale_y_continuous
#' scale_x_continuous aes_string margin geom_errorbar element_text
#' @importFrom graphics arrows axis lines mtext plot legend
#' @export
#' @examples
#' data(PEAC_minimal_load)
#' 
#' disp <- apply(tpm, 1, function(x){ 
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2) 
#' })
#' 
#' MS4A1Fit <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      id = 'PATID',
#'                      countdata = tpm['MS4A1', ],
#'                      metadata = metadata,
#'                      dispersion = disp['MS4A1'],
#'                      verbose=FALSE)
#'                      
#' modelPlot(MS4A1Fit,
#'           "MS4A1", 
#'           x1Label="Timepoint", 
#'           x2Label="EULAR_6m", 
#'           colours = c('skyblue', 'goldenrod1', 'mediumvioletred'),
#'           xTitle="Time")


modelPlot <- function(glmmResult,
                      geneName,
                      x1Label="Timepoint",
                      x2Label,
                      xTitle=NULL,
                      yTitle="Gene Expression",
                      title=NULL,
                      logTransform=FALSE,
                      shapes=21,
                      colours=c('blue'),
                      x2Offset=6, 
                      lineWidth=1, 
                      markerSize=5,
                      fontSize=10,
                      overlap=TRUE,
                      addErrorbars=TRUE, 
                      graphics="ggplot", 
                      ...) {
  
  if(! graphics %in% c("ggplot", "base")){
    stop("graphics must be either 'ggplot' or 'base'")
  }
  if(! geneName %in% rownames(glmmResult@predict)){
    stop("geneName must be in rownames(glmmResult@predict)")
  }
  if(ncol(glmmResult@modelData) != 2){
    stop(paste("These plots only work for interactions between two variable.", 
               "Therefore nrow(glmmResult@modelData) should be 2."))
  }
  
  
  # Set up the plotting data
  modelData <- glmmResult@modelData
  outLabels <- apply(glmmResult@modelData, 1, 
                     function(x) paste(x, collapse="_"))
  modelData$y <- glmmResult@predict[geneName, paste0("y_",  outLabels)]
  modelData$LCI <- glmmResult@predict[geneName, paste0("LCI_",  outLabels)]
  modelData$UCI <- glmmResult@predict[geneName, paste0("UCI_",  outLabels)]
  x2Values <- levels(modelData[, x2Label])
  x1Values <- levels(factor(modelData[, x1Label]))
  modelData[, x2Label] <- as.numeric(modelData[, x2Label])
  
  # Define the plotting parameters
  maxX2 <- max(as.numeric(modelData[, x2Label]))
  maxX1 <- max(as.numeric(modelData[, x1Label]))
  if (length(shapes)==1) shapes <- rep(shapes, maxX2)
  if (length(colours)==1) colours <- rep(colours, maxX2)
  yCover <- glmmResult@predict[geneName, ]
  yLim <- range(yCover[is.finite(yCover) & yCover != 0])
  x <- as.numeric(modelData[,x1Label]) +
    (as.numeric(modelData[,x2Label])-1) * maxX1
  p <- glmmResult@stats[geneName, grepl("P_", colnames(glmmResult@stats))]
  pval <- vapply(p, format, FUN.VALUE="1", digits=2)
  modelData$group <- factor(modelData[, x2Label], 
                            labels=levels(factor(
                              glmmResult@metadata[, x2Label])))
  modelData$x1Numeric <- as.numeric(modelData[, x1Label])
  modelData$x <- x
  
  if(is.null(title)) title <- geneName
  
  if(addErrorbars){
    minLabel <- min(c(modelData$LCI, modelData$y), na.rm=TRUE)
  } else { minLabel <- min(modelData$y, na.rm=TRUE) } 
  
  # Plot base graphics
  if(graphics == "base"){
    if(logTransform) log <- "y" else log <- ""
    if(overlap) xUse <- as.numeric(factor(modelData[, x1Label])) else xUse <- x
    
    plot(xUse, modelData$y, type='p', bty='l', las=1, xaxt='n',
         cex.axis=fontSize, cex.lab=fontSize,
         cex=markerSize,
         pch=shapes[modelData[,x2Label]], 
         bg=colours[modelData[,x2Label]], 
         xlab=xTitle,
         ylim=yLim, ylab=yTitle, log=log,
         ...,
         panel.first={
           for (i in 1:maxX2) {
             lines(xUse[modelData[,2]==i], 
                   modelData$y[modelData[,2]==i],
                   col=colours[i], lwd=2)
           }
           if(addErrorbars){
             for (i in 1:maxX2) {
               arrows(x0=xUse[modelData[,x2Label]==i],
                      y0=modelData$LCI[modelData[,x2Label]==i],
                      y1=modelData$UCI[modelData[,x2Label]==i],
                      angle=90, code=3, length=0.08, col=colours[i])
             }}
         })
    
    if (length(x2Values)==3) x2Values <- gsub('\\..*', '', x2Values)
    if(! overlap){
      axis(1, length(x1Values)*(seq_along(x2Values)-1) +
             length(x1Values)/2+0.5, 
           labels=levels(factor(glmmResult@metadata[, x2Label])), 
           line=1.5, cex.axis=fontSize, tick=FALSE)
      axis(1, modelData$x, labels=modelData[, x1Label], cex.axis=fontSize)
    } else{
      axis(1, unique(xUse), labels=levels(factor(modelData[, x1Label])), 
           cex.axis=fontSize)
      legend("top", legend=levels(modelData$group),
             col=colours, lty=1, cex=fontSize, box.lwd=0)
    }
    if(title!="") mtext(title, side=3, adj=0, padj=-3, cex=fontSize)
    
    mtext(bquote(
      paste("P"[.(x1Label)]*"=", .(pval[1]),
            ", P"[.(x2Label)]*"=", .(pval[2]),
            ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
            .(pval[3]))),
      side=3, adj=0.04, cex=fontSize)
    
    # Plot ggplot
  } else{
    
    if(overlap){
      
      p <- ggplot(modelData, aes_string(x="x1Numeric", y="y", fill="group", 
                                        shape="group", group="group", 
                                        color="group")) +
        geom_line(size=lineWidth) +
        geom_point(size=markerSize) +
        scale_x_continuous(labels=levels(modelData[, x1Label]), 
                           breaks=unique(modelData$x1Numeric)) +
        theme_classic() +
        theme(plot.margin=margin(7,0,14,0), text=element_text(size=fontSize), 
              ...) 
    } else{
      p <- ggplot(modelData, aes_string(x="x", y="y", fill="group", 
                                        shape="group", color="group")) +
        geom_line(size=lineWidth) +
        geom_point(size=markerSize) +
        scale_x_continuous(labels=modelData[, x1Label], 
                           breaks=modelData$x, name="") +
        geom_text(data=data.frame(
          label=levels(factor(glmmResult@metadata[, x2Label])),
          x=length(x1Values)*(seq_along(x2Values)-1) +
            length(x1Values)/2+0.5,
          y = minLabel), size=rel(3),
          mapping=aes_string(label="label", x="x", y="y"), hjust = 0.5,
          vjust=x2Offset, inherit.aes=FALSE)  +
        
        theme_classic() +
        theme(legend.position = "none", plot.margin=margin(7,0,14,0), 
              text=element_text(size=fontSize), ...)  +
        coord_cartesian(clip = 'off', expand=TRUE)
    }
    p <- p +
      scale_fill_manual(values=colours, name=x2Label) +
      scale_color_manual(values=colours, name=x2Label) +
      scale_shape_manual(values=shapes, name=x2Label) +
      labs(subtitle= bquote(
        paste("P"[.(x1Label)]*"=", .(pval[1]),
              ", P"[.(x2Label)]*"=", .(pval[2]),
              ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
              .(pval[3]))), y=yTitle, x=xTitle, title=title) 
    
    if(addErrorbars) p <- p + geom_errorbar(aes_string(ymin="LCI", ymax="UCI"), 
                                            width=.25, size=lineWidth) 
    if(logTransform) p <- p + scale_y_continuous(trans='log10')
    
    return(p)
  }
}
