#' Paired plots
#'
#' Paired plots to show differences between groups and over time
#' @param glmmResult A glmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}}
#' @param geneName The gene/row name to be plotted
#' @param x1Label The name of the first (inner) x parameter. This must be able
#' to be paired using the ID.
#' @param x2Label The name of the second (outer) x parameter
#' @param IDColumn Column name of sample IDs for pairing
#' @param xTitle Title for the x axis
#' @param yTitle Title for the y axis
#' @param title Plot title. If NULL gene name is used
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shapes (default=21)
#' @param colours The marker colours (default='red')
#' @param lineColour The line colours (default='grey60')
#' @param markerSize Size of markers (default=2)
#' @param fontSize Plot font size
#' @param alpha Line and marker opacity (default=0.7)
#' @param x2Offset Vertical adjustment to secondary x-axis (default=6)
#' @param pairedOnly Logical whether to only plot paired samples (default=TRUE)
#' @param graphics Which graphic system to use (options = "base" or "ggplot")
#' @param addModel Whether to add the fit model with markers (default=TRUE)
#' @param modelSize Size of model points (default=3)
#' @param modelColour Colour of model fit markers (default="black")
#' @param modelLineSize Size of model points (default=1)
#' @param modelLineColour Colour of model fit lines (default="black")
#' @param addBox Logical whether to add boxplots for mean and IQR.
#' @param addViolins Logical whether to add half violin-plots (ggplot only),
#' default=TRUE
#' @param violinWidth Width of violin plots (default=0.5)
#' @param ... Other parameters to pass to
#' \code{\link[graphics:plot]{graphics::plot()}} or
#' \code{\link[ggplot2:theme]{ggplot2::theme()}}.
#' @return Returns a paired plot for matched samples.
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme geom_boxplot scale_y_continuous
#' scale_x_continuous aes_string annotate margin element_text rel
#' @importFrom graphics arrows axis lines mtext plot segments points boxplot
#' @importFrom gghalves geom_half_violin
#' @export
#' @examples
#' data(PEAC_minimal_load)
#'
#' disp <- apply(tpm, 1, function(x){
#' (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
#' })
#'
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      id = 'PATID',
#'                      countdata = tpm['MS4A1', ],
#'                      metadata = metadata,
#'                      dispersion = disp['MS4A1'],
#'                      removeDuplicatedMeasures=TRUE,
#'                      verbose=FALSE)
#'
#' pairedPlot(glmmResult=MS4A1glmm,
#'            geneName = 'MS4A1',
#'            x1Label = 'Timepoint',
#'            x2Label='EULAR_6m',
#'            IDColumn = 'PATID',
#'            colours = c('skyblue', 'goldenrod1', 'mediumvioletred'),
#'            graphics = 'base')

pairedPlot <- function(glmmResult,
                       geneName = NULL,
                       x1Label = NULL,
                       x2Label = NULL,
                       IDColumn = 'ID',
                       xTitle = NULL,
                       yTitle = "Gene Expression",
                       title=NULL,
                       logTransform=FALSE,
                       shapes=21,
                       colours='red',
                       lineColour='grey60',
                       markerSize=2,
                       fontSize=NULL,
                       alpha=0.7,
                       x2Offset=6,
                       pairedOnly=TRUE,
                       graphics="base",
                       addModel=TRUE,
                       modelSize=3,
                       modelLineSize=1,
                       modelColour="black",
                       addBox=FALSE,
                       addViolins=TRUE,
                       violinWidth=0.5,
                       ...) {

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
    if(! geneName %in% rownames(glmmResult@countdata)){
      stop("geneName must be in rownames(glmmResult@countdata)")
    }
    if(ncol(glmmResult@modelData) != 2){
      stop(
        paste("These plots only work for interactions between two variable.",
              "Therefore nrow(glmmResult@modelData) should be 2."))
    }

    # Set up plotting data frame
    df_long <- cbind(
      glmmResult@metadata[, c(IDColumn, x1Label, x2Label)],
      geneExp=as.numeric(glmmResult@countdata[geneName, ]))
    pval <- glmmResult@stats[geneName, grepl("P_", colnames(glmmResult@stats))]
    pval <- vapply(pval, format, FUN.VALUE="character", digits=2)
    x1Values <- levels(droplevels(factor(glmmResult@modelData[, x1Label])))
    x2Values <- levels(droplevels(factor(glmmResult@modelData[, x2Label])))
    labs <- levels(droplevels(factor(glmmResult@modelData[, x2Label])))
    df_long$x2 <- as.numeric(factor(df_long[, x2Label]))
    df_long$x1 <- as.numeric(factor(df_long[, x1Label]))
    maxX2 <- max(as.numeric(df_long$x2), na.rm=TRUE)
    df_long$x1Factors <- as.numeric(factor(df_long[, x1Label])) +
      (df_long$x2-1) * length(unique(df_long[, x1Label]))
    df_long$id <- df_long[, glmmResult@variables]


    # Set up model fit data
    modelData <- glmmResult@modelData
    df_mean <- glmmResult@predict[geneName, ]
    outLabels <- apply(modelData, 1, function(x) paste(x, collapse="_"))
    df_mean <- data.frame("x1Factors"=seq_along(modelData[, 1]),
                          "y"= df_mean[ paste0("y_",outLabels)],
                          "lower" =  df_mean[paste0("LCI_", outLabels)],
                          "upper" =  df_mean[paste0("UCI_", outLabels)],
                          "group" = glmmResult@modelData[, x2Label])

  } else{
    stop("glmmResult must be an output from either glmmGene or glmmSeq")
  }

  # Check for duplicated measurements
  checkColumns <-  c(IDColumn, x1Label, x2Label)
  duplicatedMeasures <- df_long[duplicated(df_long[, checkColumns]), ]

  if(nrow(duplicatedMeasures) > 0){
    keepMeasures <- apply(df_long[, checkColumns], 1, function(x) {
      ! paste(x, collapse="-") %in%
        apply(duplicatedMeasures[ , checkColumns], 1 , paste , collapse = "-" )
    })
    df_long <- df_long[keepMeasures, ]
    str <- apply(duplicatedMeasures[, c(x1Label, x2Label)], 1, function(x){
      paste0(x1Label, "=", x[1], "; ", x2Label, "=", x[2])
    })

    warning(paste0(paste(duplicatedMeasures[, IDColumn], collapse=", "),
                   " has multiple entries for: ", str,
                   ". These will all be removed from plotting."))
  }

  # Convert to wide format
  df <- stats::reshape(df_long[, c('id', 'x1', 'geneExp')], timevar='x1',
                       idvar='id', v.names='geneExp', direction='wide')


  if (length(shapes)==1) shapes <- rep(shapes, maxX2)
  if (length(colours)==1) colours <- rep(colours, maxX2)

  # keep only paired samples
  paired <- df$id[complete.cases(df)]
  if (pairedOnly) df_long <- df_long[df_long$id %in% paired, ]

  df_long <- df_long[complete.cases(df_long), ]

  if(is.null(title)) title <- geneName

  # Generate base plots
  if(graphics == "base"){
    if(logTransform) log <- "y" else log <- ""
    if(is.null(xTitle)) xTitle <- NA

    plot(as.numeric(df_long$x1Factors), df_long$geneExp,
         type='p', bty='l', las=2,
         xaxt='n', cex.axis=fontSize, cex.lab=fontSize,
         pch=shapes[df_long$x2], bg=colours[df_long$x2], col = colours[df_long$x2],
         cex=markerSize, xlab=xTitle, ylab=yTitle,
         log=log,
         ...,
         panel.first={
           for (i in unique(df_long$id)) {
             lines(df_long$x1Factors[df_long$id==i],
                   df_long$geneExp[df_long$id==i],
                   col=lineColour[df_long[which(df_long$id == i), x2Label]]))}
         })

    if(addModel){
      for(i in unique(df_mean$group)){
        lines(df_mean$x1Factors[df_mean$group == i],
              df_mean$y[df_mean$group == i],
              lwd=modelSize+1, col=modelColour[i])
      }
      segments(df_mean$x1Factors, df_mean$upper,
               df_mean$x1Factors, df_mean$lower,
               lwd=modelSize+1, col=rep(modelColour, each=length(unique(df_mean$group))))
      points(df_mean$x1Factors, df_mean$y, type = "p", bg=rep(modelColour, each=length(unique(df_mean$group))),
             col=rep(modelColour, each=length(unique(df_mean$group))), pch=21, cex=modelSize)
    }

    if(addBox){
      boxplot(geneExp ~ x1Factors , data = df_long, outcex=0, boxwex=0.25,
              add=TRUE, col=NA, frame = FALSE, axes=FALSE)
    }
    axis(1, seq_along(modelData[, 1]), labels=rep(x1Values, length(x2Values)),
         cex.axis=fontSize)
    axis(1, length(x1Values)*(seq_along(x2Values)-1) +
           length(x1Values)/2+0.5, labels=labs,
         line=1.5, cex.axis=fontSize, tick=FALSE)
    if(title!="") mtext(title, side=3, adj=0, padj=-3, cex=fontSize)

    mtext(bquote(
      paste("P"[.(x1Label)]*"=", .(pval[1]),
            ", P"[.(x2Label)]*"=", .(pval[2]),
            ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
            .(pval[3]))), cex=fontSize,
      side=3, adj=0)


  } else{ # Generate ggplot plots
    df_long$x1Factors <- factor(df_long$x1Factors)
    df_long$x2 <- factor(df_long$x2)

    p <- ggplot(df_long, aes_string(x="x1Factors", y="geneExp", group="id",
                                    shape="x2", fill="x2", color="x2")) +
      geom_line(color=lineColour, alpha=alpha) +
      geom_point(size=markerSize, alpha=alpha) +
      theme_classic() +
      scale_fill_manual(values=colours) +
      scale_color_manual(values=colours) +
      scale_shape_manual(values=shapes) +
      labs(x=xTitle, y=yTitle, title=title) +
      theme(legend.position = "none", plot.margin=margin(7,0,14+x2Offset,0),
            text=element_text(size=fontSize), ...) +
      scale_x_discrete(labels=setNames(rep(x1Values, length(x2Values)),
                                       seq_along(modelData[, 1]))) +
      geom_text(data=data.frame(label = labs,
                                x=length(x1Values)*(seq_along(x2Values)-1) +
                                  length(x1Values)/2+0.5,
                                y = min(c(modelData$LCI, df_long$geneExp),
                                        na.rm=TRUE)),
                size=rel(3),
                mapping=aes_string(label="label", x="x", y="y"), hjust = 0.5,
                nudge_x=0, vjust=x2Offset, inherit.aes=FALSE)  +
      coord_cartesian(clip = 'off')

    if(addBox) {
      p <- p + geom_boxplot(mapping=aes_string(x="x1Factors", y="geneExp"),
                            inherit.aes=FALSE,
                            alpha=alpha*0.7, outlier.shape=NA, width=0.2)
    }

    if(addViolins){
      df_long$nudged <- as.numeric(df_long$x1Factors) +
        0.06*ifelse(df_long$x1 %% 2 == 0, 1, -1)
      df_left <- df_long[df_long$x1 %% 2 != 0, ]
      lLength <- length(which(table(df_left$nudged)>2))
      leftColours <- rep(
        colours,
        each=length(unique(df_left$x1)))[table(df_left$nudged)>2]
      df_right <- df_long[df_long$x1 %% 2 == 0, ]
      rLength <- length(which(table(df_right$nudged)>2))
      rightColours <- rep(
        colours,
        each=length(unique(df_right$x1)))[table(df_right$nudged)>2]

      p <- p +
        # Violins on the right
        geom_half_violin(data = df_right,
                         inherit.aes = FALSE, alpha=alpha,
                         mapping=aes_string(x="nudged",
                                            y="geneExp", group="x1Factors"),
                         side=c("r"), width=violinWidth,
                         color=rep(rightColours,
                                   each=512*rLength/length(rightColours)),
                         fill=rep(rightColours,
                                  each=512*rLength/length(rightColours)))  +
        # Violins on the left
        geom_half_violin(data = df_left,
                         inherit.aes = FALSE, alpha=alpha,
                         mapping=aes_string(x="nudged",
                                            y="geneExp", group="x1Factors"),
                         side=c("l"), width=violinWidth,
                         color=rep(leftColours,
                                   each=512*lLength/length(leftColours)),
                         fill=rep(leftColours,
                                  each=512*lLength/length(leftColours)))
    }

    if(addModel & class(glmmResult) == "GlmmSeq"){
      p <- p +
        annotate("line", x = df_mean$x1Factors, y = df_mean$y,
                 group=df_mean$group, size=modelLineSize,
                 color=modelLineColour) +
        annotate("errorbar", x = df_mean$x1Factors, y = df_mean$y,
                 color=modelLineColour, ymin=df_mean$lower, ymax=df_mean$upper,
                 width=0.2, size=modelLineSize) +
        annotate("point", x = df_mean$x1Factors, y = df_mean$y,
                 size=modelSize, color=modelColour)
    }


    if(logTransform) p <- p + scale_y_continuous(trans='log10')


    p <- p + labs(subtitle= bquote(
      paste("P"[.(x1Label)]*"=", .(pval[1]),
            ", P"[.(x2Label)]*"=", .(pval[2]),
            ", P"[paste(.(x1Label), ":", .(x2Label))]*"=",
            .(pval[3]))))

    return(p)
  }

}
