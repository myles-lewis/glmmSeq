#' Model plot
#'
#' Model plots to show the overall differences between groups and over time
#' @param result A glmmSeq object
#' @param gene Gene name to plot
#' @param pch The marker shape (default=21)
#' @param bg The marker colour, default=c('red', 'blue')
#' @param graphics Which graphic system to use (options = "base" or "ggplot")
#' @return Returns a model plot
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme scale_y_continuous
#' scale_x_continuous aes_string
#' @importFrom graphics arrows axis lines mtext plot
#' @keywords hplot
#' @export

modelPlot <- function(result,
                      gene,
                      pch=21,
                      bg=c('red', 'blue'),
                      graphics="ggplot") {

  if(! graphics %in% c("ggplot", "base")){
    stop("graphics must be either 'ggplot' or 'base'")
  }
  if(! gene %in% rownames(result@predict)){
    stop("gene must be in rownames(result@predict)")
  }

  modelData <- result@modelData
  s <- nrow(modelData)
  modelData$y <- result@predict[gene, 1:s]
  modelData$LCI <- result@predict[gene, 1:s + s]
  modelData$UCI <- result@predict[gene, 1:s + s*2]
  labs <- levels(modelData[,2])
  modelData[,2] <- as.numeric(modelData[,2])
  maxY <- max(as.numeric(modelData[,2]))
  maxX <- max(as.numeric(modelData[,1]))
  if (length(pch)==1) pch <- rep(pch, maxY)
  if (length(bg)==1) bg <- rep(bg, maxY)
  yCover <- result@predict[gene, ]
  yLim <- range(yCover[is.finite(yCover) & yCover != 0])
  x <- as.numeric(modelData[,1]) + (as.numeric(modelData[,2])-1) * maxX
  p <- result@stats[gene, grepl("P_", colnames(result@stats))]
  pval <- vapply(p, format, FUN.VALUE=1, digits=2)

  if(plot == "base"){
    plot(x, modelData$y, type='p', bty='l', las=1, xaxt='n',
         cex.axis=1.5, cex.lab=1.8,
         pch=pch[modelData[,2]], bg=bg[modelData[,2]], cex=1.8, xlab=NA,
         ylim=yLim, ylab=gene, log="y",
         panel.first={
           for (i in 1:maxY) {
             lines(x[modelData[,2]==i], modelData$y[modelData[,2]==i],
                   col=bg[i], lwd=2)
           }
           arrows(x0=x, y0=modelData$LCI, y1=modelData$UCI,
                  angle=90, code=3, length=0.08)
         })
    axis(1, 1:(2*maxY), labels=rep(modelData$Timepoint, maxY), cex.axis=1.6)
    if (length(labs)==3) labs <- gsub('\\..*', '', labs)
    axis(1, 1:maxY*2-0.5, labels=labs, line=1.5, cex.axis=1.4, tick=FALSE)

    for (i in 1:seq_along(pval)) {
      mtext(bquote( .(paste0(names(pval[i]), " = ", pval[i])) ),
            line=length(pval)-i, side=3, adj=0, cex=1)
    }
  } else{
    df <- data.frame(x=x, y=modelData$y, c=factor(modelData[, 2]))
    ggplot(df, aes_string(x="x", y="y", fill="c", shape="c"), colour="black") +
      geom_point(size=3) +
      geom_line(color="black") +
      scale_fill_manual(values=bg) +
      scale_shape_manual(values=pch) +

      scale_x_continuous(labels=setNames(modelData$Timepoint, 1:4)) +
      labs(title= bquote(.(paste(paste(names(pval), "=", pval),
                                 collapse="\n"))),
           y=gene, x="") +
      theme_classic() +
      theme(legend.position = "none") +
      geom_text(data=data.frame(label = unique(modelData$erosionstatus),
                                x = c(1.5, 3.5),
                                y = 10),
                mapping=aes_string(label="label", x="x", y="y"), hjust = 0.5,
                vjust=5, nudge_x=0,  inherit.aes=FALSE)  +



      coord_cartesian(clip = 'off') +
      scale_y_continuous(trans='log10')
  }
}
