#' Paired plots
#'
#' Paired plots to show differences between groups and over time
#' @param glmmResult A GlmmSeq object from \code{glmmSeq}.
#' @param geneName The gene name to be plotted.
#' @param x1Label The name of the first (inner) x parameter. This must be able 
#' to be paired using the ID. 
#' @param x2Label The name of the second (outer) x parameter
#' @param xTitle Title for the x axis
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shapes (default=21)
#' @param colours The marker colours (default='red')
#' @param lineColour The line colours (default='grey60')
#' @param alpha Line and marker opacity (default=0.7)
#' @param pairedOnly Logical whether to only plot paired samples (default=TRUE)
#' @param title Plot title
#' @param graphics Which graphic system to use (options = "base" or "ggplot")
#' @param addModel Whether to add the fit model with markers
#' @param modSize Size of model points and line
#' @param modColour Colour of model fit markers
#' @param modLineColour Colour of model fit lines
#' @param addBox Logical whether to add boxplots.
#' @param addViolins Logical whether to add half violins (ggplot only)
#' @param violinWidth Width of violin plots
#' @return Returns a paired plot
#' @importFrom ggpubr compare_means ggboxplot stat_pvalue_manual
#' stat_compare_means
#' @importFrom lme4 lmer isSingular glmerControl glmer
#' @importFrom MASS negative.binomial
#' @importFrom car Anova
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme geom_boxplot scale_y_continuous
#' scale_x_continuous aes_string annotate margin
#' @importFrom graphics arrows axis lines mtext plot segments points boxplot
#' @importFrom gghalves geom_half_violin
#' @export

pairedPlot <- function(glmmResult,
                       geneName = NULL,
                       x1Label = NULL,
                       x2Label = NULL,
                       IDColumn = 'ID',
                       xTitle = NULL,
                       logTransform=FALSE,
                       shapes=21,
                       colours='red',
                       lineColour='grey60',
                       alpha=0.7,
                       pairedOnly=TRUE,
                       title="",
                       graphics="base",
                       addModel=TRUE,
                       modSize=1.5,
                       modColour="black",
                       modLineColour="black",
                       addBox=FALSE,
                       addViolins=TRUE, 
                       violinWidth=0.5,
                       ...) {
  
  suffix <- ""
  if(class(alpha) != "numeric" | alpha > 1 | alpha < 0) {
    stop("alpha must be a numer between 1 and 0")
  }
  
  # For outputs from glmmGene
  if(class(glmmResult) == "glmerMod"){
    if(! x1Label %in% colnames(glmmResult@frame)) {
      stop("x1Label must be a column name in glmmResult@frame")
    }
    if(! x2Label %in% colnames(glmmResult@frame)) {
      stop("x2Label must be a column name in glmmResult@frame")
    }
    if(! IDColumn %in% colnames(glmmResult@frame)) {
      stop("IDColumn must be a column name in glmmResult@frame")
    }
    
    x1Values <-  levels(droplevels(factor(glmmResult@frame[, x1Label])))
    labs <- levels(droplevels(factor(glmmResult@frame[, x2Label])))
    
    # For outputs from glmmSeq
  } else if(class(glmmResult) == "GlmmSeq"){
    if(! x1Label %in% colnames(glmmResult@modelData)) {
      stop("x1Label must be a column name in glmmResult@modelData")
    }
    if(! x2Label %in% colnames(glmmResult@modelData)) {
      stop("x2Label must be a column name in glmmResult@modelData")
    }
    if(! IDColumn %in% colnames(glmmResult@metadata)) {
      stop("IDColumn must be a column name in glmmResult@metadata")
    }
    
    # Set up plotting data frame
    df_long <- cbind(
      glmmResult@metadata[, c(IDColumn, x1Label, x2Label)],
      gene_exp=as.numeric(glmmResult@countdata[geneName, ]))
    pval <- glmmResult@stats[geneName, grepl("P_", colnames(glmmResult@stats))]
    pval <- vapply(pval, format, FUN.VALUE="character", digits=2)
    x1Values <- levels(droplevels(factor(glmmResult@modelData[, x1Label])))
    labs <- levels(droplevels(factor(glmmResult@modelData[, x2Label])))
    df_long$x2 <- as.numeric(factor(df_long[, x2Label]))
    df_long$x1 <- as.numeric(factor(df_long[, x1Label]))
    maxX2 <- max(as.numeric(df_long$x2), na.rm=TRUE)
    df_long$x1Factors <- as.numeric(factor(df_long[, x1Label])) + 
      (df_long$x2-1) * length(unique(df_long[, x1Label]))
    df_long$geneExp <- as.numeric(glmmResult@countdata[geneName, ])
    df_long$id <- df_long[, glmmResult@variables]
    
    # Catch if singular
    if (glmmResult@optInfo[geneName, "Singular"] == 1) suffix <- "s"
    
    # Set up model fit data
    df_mean <- glmmResult@predict[geneName, ]
    df_mean <- data.frame("x1Factors"=seq_along(glmmResult@modelData[, 1]),
                          "y"= df_mean[
                            paste0("y",
                                   seq_along(glmmResult@modelData[, 1]))],
                          "lower" = 
                            df_mean[
                              paste0("LCI", 
                                     seq_along(glmmResult@modelData[, 1]))],
                          "upper" = 
                            df_mean[
                              paste0("UCI", 
                                     seq_along(glmmResult@modelData[, 1]))],
                          "id" = rep(seq_along(unique(df_long$x2)),
                                     each=length(unique(df_long$x2))))
    
  } else{
    stop("glmmResult must be an output from either glmmGene or glmmSeq")
  }
  
  # Convert to wide format
  df <- reshape(df_long[, c('id', 'x1', 'geneExp')], timevar='x1',
                idvar='id', v.names='geneExp', direction='wide')
  
  
  if (length(shapes)==1) shapes <- rep(shapes, maxX2)
  if (length(colours)==1) colours <- rep(colours, maxX2)
  
  # keep only paired samples
  paired <- df$id[complete.cases(df)] 
  if (pairedOnly) df_long <- df_long[df_long$id %in% paired, ]
  
  df_long <- df_long[complete.cases(df_long), ]
  
  
  # Generate base plots
  if(graphics == "base"){
    if(logTransform) log <- "y" else log <- ""
    if(is.null(xTitle)) xTitle <- NA
    
    plot(as.numeric(df_long$x1Factors), df_long$geneExp +.1,
         type='p', bty='l', las=1,
         xaxt='n', cex.axis=1.5, cex.lab=1.8,
         pch=shapes[df_long$x2], bg=colours[df_long$x2],
         cex=1.8, xlab=xTitle, ylab=geneName,
         log=log, 
         ...,
         panel.first={
           for (i in paired) {lines(df_long$x1Factors[df_long$id==i],
                                    df_long$geneExp[df_long$id==i]+.1, 
                                    col=lineColour)}
         })
    
    if(addModel){
      for(i in unique(df_mean$id)){
        lines(df_mean$x1Factors[df_mean$id == i], df_mean$y[df_mean$id == i], 
              lwd=modSize+1, col=modLineColour)
      }
      segments(df_mean$x1Factors, df_mean$upper, 
               df_mean$x1Factors, df_mean$lower, 
               lwd=modSize+1, col=modLineColour)
      points(df_mean$x1Factors, df_mean$y, type = "p", bg=modColour, 
             col=modLineColour, pch=21, cex=modSize)
    }
    
    if(addBox){
      boxplot(geneExp + 0.1 ~ x1Factors , data = df_long, 
              add=TRUE, col=NA, frame = FALSE, axes=FALSE)
    }
    axis(1, 1:(2*maxX2), labels=rep(x1Values, maxX2), cex.axis=1.6)
    if (length(labs)==3) labs <- gsub('\\..*', '', labs)
    axis(1, 1:maxX2*2-0.5, labels=labs, line=1.5, cex.axis=1.4, tick=FALSE)
    if(title!="") mtext(title, side=3, adj=0, padj=-3, cex=1.5)
    if (suffix!="") suffix <- paste0(" [", suffix, "]")
    
    mtext(bquote(
      paste("P"[.(x1Label)]*"=", .(pval[1]),
            ", P"[.(x2Label)]*"=", .(pval[2]),
            ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
            .(pval[3]), .(suffix))), 
      side=3, adj=0.04)
    
    
  } else{ # Generate ggplot plots
    df_long$x1Factors <- factor(df_long$x1Factors)
    df_long$x2 <- factor(df_long$x2)
    
    p <- ggplot(df_long, aes_string(x="x1Factors", y="geneExp", group="id",
                                    shape="x2", fill="x2", color="x2")) +
      geom_line(color=lineColour, alpha=alpha) +
      geom_point(size=3, alpha=alpha) +
      theme_classic() +
      scale_fill_manual(values=colours) +
      scale_color_manual(values=colours) +
      scale_shape_manual(values=shapes) +
      labs(x=xTitle, y=geneName, title=title) +
      theme(legend.position = "none", plot.margin=margin(0,0,15,0)) +
      scale_x_discrete(labels=setNames(rep(x1Values, 2), 1:4)) +
      geom_text(data=data.frame(label = labs, x = c(1.5, 3.5),
                                y = min(df_long$geneExp)),
                mapping=aes_string(label="label", x="x", y="y"), hjust = 0.5,
                vjust=5, nudge_x=0,  inherit.aes=FALSE)  +
      coord_cartesian(clip = 'off')
    
    if(addBox) {
      p <- p + geom_boxplot(mapping=aes_string(x="x1Factors", y="y"), 
                            inherit.aes=FALSE,
                            alpha=alpha*0.7, outlier.shape=NA, width=0.2, 
                            fill=rep(colours, each=2))
    }
    
    if(addViolins){
      df_long$nudged <- as.numeric(df_long$x1Factors) + 
        0.06*ifelse(df_long$x1 %% 2 == 0, 1, -1)
      
      p <- p +
        # Violins on the right
        geom_half_violin(data = df_long[df_long$x1 %% 2 == 0, ], 
                         inherit.aes = FALSE, alpha=alpha,
                         mapping=aes_string(x="nudged", 
                                            y="geneExp", group="x1Factors"), 
                         side=c("r"), color=rep(colours, each=1024/2), 
                         fill=rep(colours, each=1024/2), width=violinWidth)  + 
        # Violins on the left 
        geom_half_violin(data = df_long[df_long$x1 %% 2 != 0, ],
                         inherit.aes = FALSE, alpha=alpha,
                         mapping=aes_string(x="nudged", 
                                            y="geneExp", group="x1Factors"), 
                         side=c("l"), color=rep(colours, each=1024/2), 
                         fill=rep(colours, each=1024/2), width=violinWidth) 
    }
    
    if(addModel & class(glmmResult) == "GlmmSeq"){
      p <- p +
        annotate("line", x = df_mean$x1Factors, y = df_mean$y,
                 group=df_mean$id, size=modSize, color=modLineColour) +
        annotate("errorbar", x = df_mean$x1Factors, y = df_mean$y,
                 color=modLineColour, ymin=df_mean$lower, ymax=df_mean$upper,
                 width=0.2, size=modSize) +
        annotate("point", x = df_mean$x1Factors, y = df_mean$y,
                 size=modSize+3, color=modColour)
    }
    
    
    if(logTransform) p <- p + scale_y_continuous(trans='log10')
    
    if (suffix!="") suffix <- paste0(" [", suffix, "]")
    
    p <- p + labs(subtitle= bquote(
      paste("P"[.(x1Label)]*"=", .(pval[1]),
            ", P"[.(x2Label)]*"=", .(pval[2]),
            ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
            .(pval[3]), .(suffix))))
    
    return(p)  
  }
  
}
