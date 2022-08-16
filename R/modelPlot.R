#' Mixed model effects plot
#'
#' Plot to show differences between groups over time using base graphics.
#' 
#' @param object A glmmSeq/lmmSeq object created by
#' \code{\link[glmmSeq:glmmSeq]{glmmSeq::glmmSeq()}} or 
#' \code{\link[glmmSeq:lmmSeq]{glmmSeq::lmmSeq()}}
#' @param geneName The gene/row name to be plotted
#' @param x1var The name of the first (inner) x parameter, typically 'time'.
#'   This is anticipated to have different values when matched by ID.
#' @param x2var The name of an optional second (outer) x parameter, which should
#'   be a factor.
#' @param x2shift Amount to shift along x axis for each level of `x2var`. By
#'   default the function will arrange each level of `x2var` side by side. Lower
#'   values of `x2shift` or `x2shift = 0` can be used to overlap plots similar
#'   to 'dodge' or stagger them.
#' @param xlab Title for the x axis
#' @param ylab Title for the y axis
#' @param plab Optional character vector of labels for p-values. These must
#'   align with column names in `object@stats$pvals`.
#' @param title Plot title. If NULL gene name is used
#' @param logTransform Whether to perform a log10 transform on the y axis
#' @param shapes The marker shapes (default=19)
#' @param colours The marker colours (default='red') as vector or named vector
#' @param lineColours The line colours (default='grey60') as vector or named 
#' vector
#' @param markerSize Size of markers (default=2)
#' @param fontSize Plot font size
#' @param alpha Line and marker opacity (default=0.7)
#' @param addModel Whether to add the fit model with markers (default=TRUE)
#' @param addPoints Whether to add underlying data points (default=TRUE)
#' @param modelSize Size of model points (default=2)
#' @param modelColours Colour of model fit markers (default="black") as vector 
#' or named vector
#' @param modelLineSize Size of model points (default=1) as vector or named 
#' vector
#' @param modelLineColours Colour of model fit lines.
#' @param errorBarLwd Line width of error bars
#' @param errorBarLength Head width of error bars
#' @param ... Other parameters to pass to
#' \code{\link[graphics:plot]{graphics::plot()}}
#' @return Returns a paired plot for matched samples
#' @importFrom graphics arrows axis lines mtext plot segments points boxplot
#' @export
#' @examples
#' data(PEAC_minimal_load)
#'
#' disp <- apply(tpm, 1, function(x){
#'   (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
#' })
#'
#' MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
#'                      countdata = tpm[1:2, ],
#'                      metadata = metadata,
#'                      dispersion = disp)
#'
#' modelPlot(object=MS4A1glmm,
#'           geneName = 'MS4A1',
#'           x1var = 'Timepoint',
#'           x2var='EULAR_6m')

modelPlot <- function(object,
                      geneName = NULL,
                      x1var = NULL,
                      x2var = NULL,
                      x2shift = NULL,
                      xlab = NA,
                      ylab = geneName,
                      plab = NULL,
                      title = geneName,
                      logTransform = is(object, "GlmmSeq"),
                      shapes = 21,
                      colours = 'grey60',
                      lineColours = 'grey60',
                      markerSize = 0.5,
                      fontSize = NULL,
                      alpha = 0.7,
                      addModel = TRUE,
                      addPoints = TRUE,
                      modelSize = 2,
                      modelColours = "royalblue",
                      modelLineSize = 1,
                      modelLineColours = modelColours,
                      errorBarLwd = 2.5,
                      errorBarLength = 0.05,
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
  
  # Generate base plots
  log <- if(logTransform) "y" else ""
  if(addModel) {
    myYlim <- if (addPoints) {
      range(c(df_model[, c('lower', 'upper')], df_long$y))
    } else range(df_model[, c('lower', 'upper')])
  } else myYlim <- NULL
  if (addPoints) {
    plot(df_long$x, df_long$y,
         ylim = myYlim, type='p', bty='l', las=2,
         xaxt='n', cex.axis=fontSize, cex.lab=fontSize,
         pch=19, 
         col=colours[df_long$x2],
         cex=markerSize, xlab=xlab, ylab=ylab,
         log=log, xlim = xlim,
         ...,
         panel.first={
           for (i in as.character(unique(df_long$id))) {
             lines(df_long$x[df_long$id==i],
                   df_long$y[df_long$id==i],
                   col=lineColours[df_long$x2[df_long$id == i]] )}
         })
  } else {
    plot(as.numeric(df_long$x), df_long$y,
         ylim = myYlim, type='n', bty='l', las=2,
         xaxt='n', cex.axis=fontSize, cex.lab=fontSize,
         xlab=xlab, ylab=ylab,
         log=log,
         xlim = xlim,
         ...)
  }
  if(addModel){
    for(i in 1:nlevels(df_model$group)){
      ind <- as.numeric(df_model$group) == i
      lines(df_model$x[ind],
            df_model$y[ind],
            lwd=modelSize+1, col=modelLineColours[i])
      arrows(df_model$x[ind], 
             df_model$upper[ind],
             df_model$x[ind], 
             df_model$lower[ind],
             lwd=errorBarLwd, col=modelLineColours[i],
             angle = 90, code = 3, length = errorBarLength)
      points(df_model$x[ind], 
             df_model$y[ind], type = "p",
             col=ifelse(shapes[i] >= 21, "black",
                        modelColours[i]),
             bg=ifelse(shapes[i] < 21, NULL,
                       modelColours[i]),
             pch=shapes[i], cex=modelSize)
    }
  }
  
  if (x2shift > xdiff) {
    axis(1, modelData[, x1var] + (as.numeric(modelData[, x2var])-1) * x2shift, 
         labels=modelData[, x1var], cex.axis=fontSize)
    axis(1, x2shift*(seq_along(x2labs)-1) + xdiff/2, labels=x2labs,
         line=1, cex.axis=fontSize, tick=FALSE)
  } else {
    axis(1, modelData[, x1var], 
         labels=modelData[, x1var], cex.axis=fontSize)
  }
  if(title!="") mtext(title, side=3, adj=0, padj=-3, cex=fontSize)
  
  if (is.null(plab)) plab <- colnames(pval)
  
  ptext <- lapply(1:ncol(pval), function(i) {
    bquote("P" [.(plab[i])] *"="* .(pval[,i]))
  })
  ptext <- bquote(.(paste(unlist(ptext), collapse = '*", "*')))
  mtext(parse(text=ptext), side=3, adj=0, cex=fontSize)
  
}


formPlot <- function(object, geneName, x1var, x2var, x2shift) {
  if (!x1var %in% colnames(object@modelData)) {
    stop("x1var must be a column name in object@modelData")}
  if (!is.null(x2var)) if (!x2var %in% colnames(object@modelData)) {
    stop("x2var must be a column name in object@modelData")}
  # if(ncol(object@modelData) > 2){
  #   stop("More than 2 variables in modelData")}
  
  maindata <- if (inherits(object, "GlmmSeq")) {
    object@countdata} else object@maindata
  if(! geneName %in% rownames(maindata)) {
    stop("geneName not found")}
  
  # Set up plotting data frame
  IDColumn <- object@vars$id
  id <- object@metadata[, IDColumn]
  y <- maindata[geneName, ]
  x <- object@metadata[, x1var]
  xdiff <- diff(range(x, na.rm = TRUE))
  if (!is.null(x2var)) {
    x2 <- as.numeric(factor(object@metadata[, x2var]))
    nsegments <- length(unique(x)) -1
    if (is.null(x2shift)) {
      x2shift <- max(x, na.rm = TRUE) + xdiff / nsegments
    }
    x <- x + (x2-1) * x2shift
  } else {
    x2 <- 1
    x2shift <- -Inf
  }
  df_long <- data.frame(id, y, x, x2)
  
  # Set up model fit data
  modelData <- object@modelData
  preds <- object@predict[geneName, ]
  s <- nrow(modelData)
  modelx <- if (!is.null(x2var)) {
    modelData[, x1var] + (as.numeric(modelData[, x2var])-1) * x2shift
  } else modelData[, x1var]
  df_model <- data.frame(x = modelx,
                         y = preds[1:s],
                         lower = preds[1:s +s],
                         upper = preds[1:s +s*2],
                         group = modelData[, x2var])
  if (is.null(x2var)) df_model$group <- 1
  
  return(list(df_long, df_model, xdiff, x2shift))
}

