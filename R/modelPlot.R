#' Model plot
#'
#' Model plots to show the overall differences between groups and over time
#' @param glmmResult A glmmSeq object
#' @param geneName Gene name to plot
#' @param x1Label The name of the first (inner) x parameter
#' @param x2Label The name of the second (outer) x parameter
#' @param xTitle Title for the x axis
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shape (default=21)
#' @param colours The marker colour, default=c('red', 'blue')
#' @param graphics Which graphic system to use (options = "base" or "ggplot")
#' @return Returns a model plot
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme scale_y_continuous
#' scale_x_continuous aes_string margin
#' @importFrom graphics arrows axis lines mtext plot
#' @export

modelPlot <- function(glmmResult,
                      geneName,
                      x1Label="Timepoint",
                      x2Label, 
                      xTitle=NULL,
                      logTransform=FALSE, 
                      shapes=21,
                      colours=c('red', 'blue'),
                      graphics="ggplot") {
  
  if(! graphics %in% c("ggplot", "base")){
    stop("graphics must be either 'ggplot' or 'base'")
  }
  if(! geneName %in% rownames(glmmResult@predict)){
    stop("geneName must be in rownames(glmmResult@predict)")
  }
  
  # Set up the plotting data
  modelData <- glmmResult@modelData
  s <- nrow(modelData)
  modelData$y <- glmmResult@predict[geneName, 1:s]
  modelData$LCI <- glmmResult@predict[geneName, 1:s + s]
  modelData$UCI <- glmmResult@predict[geneName, 1:s + s*2]
  labs <- levels(modelData[, x2Label])
  modelData[, x2Label] <- as.numeric(modelData[, x2Label])
  
  # Define the plotting parameters
  maxX2 <- max(as.numeric(modelData[, x2Label]))
  maxX1 <- max(as.numeric(modelData[, x1Label]))
  if (length(shapes)==1) shapes <- rep(shapes, maxX2)
  if (length(colours)==1) colours <- rep(colours, maxX2)
  yCover <- glmmResult@predict[geneName, ]
  yLim <- range(yCover[is.finite(yCover) & yCover != 0])
  x <- as.numeric(modelData[,x1Label]) + 
    (as.numeric(modelData[,x2Label])-1) * maxX2
  p <- glmmResult@stats[geneName, grepl("P_", colnames(glmmResult@stats))]
  pval <- vapply(p, format, FUN.VALUE="1", digits=2)
  
  # Plot base graphics
  if(graphics == "base"){
    if(logTransform) log <- "y" else log <- ""
    
    plot(x, modelData$y, type='p', bty='l', las=1, xaxt='n',
         cex.axis=1.5, cex.lab=1.8,
         pch=shapes[modelData[,2]], bg=colours[modelData[,2]], cex=1.8, 
         xlab=xTitle,
         ylim=yLim, ylab=geneName, log=log,
         panel.first={
           for (i in 1:maxX2) {
             lines(x[modelData[,2]==i], modelData$y[modelData[,2]==i],
                   col=colours[i], lwd=2)
           }
           arrows(x0=x, y0=modelData$LCI, y1=modelData$UCI,
                  angle=90, code=3, length=0.08)
         })
    axis(1, 1:(2*maxX2), labels=rep(unique(modelData[, x1Label]), maxX2), 
         cex.axis=1.6)
    if (length(labs)==3) labs <- gsub('\\..*', '', labs)
    axis(1, 1:maxX2*2-0.5, labels=labs, line=1.5, cex.axis=1.4, tick=FALSE)
    
    mtext(bquote(
      paste("P"[.(x1Label)]*"=", .(pval[1]),
            ", P"[.(x2Label)]*"=", .(pval[2]),
            ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
            .(pval[3]))), 
      side=3, adj=0.04)
    
    # Plot ggplot
  } else{
    df <- data.frame(x=x, y=modelData$y, c=factor(modelData[, 2]))
    p <- ggplot(df, aes_string(x="x", y="y", fill="c", shape="c"), 
                colour="black") +
      geom_line(color="black") +
      geom_point(size=3) +
      scale_fill_manual(values=colours) +
      scale_shape_manual(values=shapes) +
      scale_x_continuous(labels=
                           setNames(modelData[, x1Label], 
                                    seq_along(modelData[, 1])), name="") +
      labs(subtitle= bquote(
        paste("P"[.(x1Label)]*"=", .(pval[1]),
              ", P"[.(x2Label)]*"=", .(pval[2]),
              ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
              .(pval[3]))), y=geneName, x=xTitle) +
      theme_classic() +
      theme(legend.position = "none", plot.margin=margin(0,0,14,0)) +
      coord_cartesian(clip = 'off', expand=TRUE) +
      annotate("text", label = levels(factor(glmmResult@metadata[, x2Label])),
               x = c(1.5, 3.5), y = min(df$y), hjust = 0.5, vjust=5)  
    
    if(logTransform) p <- p + scale_y_continuous(trans='log10')
    
    return(p)
  }
}
