#' Paired plots
#'
#' Paired plots to show differences between groups and over time
#' @param id A character vector of subjects ids for pairing over time.
#' @param x A character vector for groups
#' @param gene_exp The gene expression
#' @param time The time point
#' @param dispersion The dispersion
#' @param datatype The type of expression data. Options include c("VST", "TPM")
#' @param plottype The plot type. Options "1"=plot all samples, "2"=plot mean
#' @param pch The marker shape (default=21)
#' @param bg='red' The marker colour (default='red')
#' @param pairedonly Logical whether to only plot paired samples (default=T)
#' @param xlab Character vector for the label on the x axis
#' @param graphics Which graphic system to use (options = "base" or "ggplot")
#' @return Returns a paired plot
#' @importFrom ggpubr compare_means ggboxplot stat_pvalue_manual
#' stat_compare_means
#' @importFrom lme4 lmer isSingular glmerControl glmer
#' @importFrom MASS negative.binomial
#' @importFrom car Anova
#' @importFrom ggplot2 ggplot geom_boxplot geom_point geom_line theme_classic
#' scale_fill_manual scale_shape_manual labs geom_text scale_x_discrete
#' coord_cartesian aes scale_color_manual theme geom_boxplot scale_y_continuous
#' @keywords hplot
#' @export

pairedPlot <- function(id,
                       time,
                       gene_exp,
                       f,
                       dispersion,
                       datatype,
                       plottype,
                       plab='drug',
                       pch=21,
                       bg='red',
                       pairedonly=T,
                       xlab="",
                       ylab="",
                       title="",
                       graphics="base",
                       addBox=TRUE,
                       ...) {

  # model is gene_exp ~ time * f + (1 | id)
  time_labs <- levels(droplevels(factor(time)))
  time <- as.numeric(factor(time))  # time
  labs <- levels(droplevels(factor(f)))
  df_long <- data.frame(id=id, time=time, gene_exp=gene_exp)

  df <- reshape(df_long[, c('id', 'time', 'gene_exp')], timevar='time',
                idvar='id', v.names='gene_exp', direction='wide')
  paired <- df$id[complete.cases(df)]  # keep only paired samples
  if (pairedonly) df_long <- df_long[df_long$id %in% paired, ]
  df_long$f <- f[match(df_long$id, id)]  # fills in missing data by ID
  df_long <- df_long[complete.cases(df_long), ]
  if (datatype=='VST') {
    fit <- lmer(gene_exp ~ time * f + (1 | id), data=df_long)
  } else if (!is.na(dispersion)) {
    fit <- try( glmer(gene_exp ~ time * f + (1 | id), data=df_long,
                      family= negative.binomial( theta=1/dispersion),
                      control=glmerControl(optimizer="bobyqa")), silent=T)
  }
  suffix <- ""
  if (class(fit) != 'try-error') {
    pval <- Anova(fit)[,3]
    pval <- sapply(pval, format, digits=2)
    if (isSingular(fit)) suffix <- "s"
    if (plottype=="2") {
      # http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
      newdata <- expand.grid(time=unique(time), f=levels(droplevels(factor(f))))
      newdata$gene_exp <- predict(fit, newdata=newdata, re.form=NA)
      designmat <- model.matrix( ~ time * f, newdata)
      predvar <- diag(designmat %*% vcov(fit) %*% t(designmat))
      newdata$SE <- sqrt(predvar)
      newdata$LCI <- newdata$gene_exp - newdata$SE * 1.96
      newdata$UCI <- newdata$gene_exp + newdata$SE * 1.96
      newdata$f <- as.numeric(newdata$f)
      newdata$time <- newdata$time + (newdata$f -1) *2
      if (datatype=='TPM') {
        newdata[, c('gene_exp', 'LCI', 'UCI')] <-
          exp(newdata[, c('gene_exp', 'LCI', 'UCI')])
      }
      gene_expcover <- unlist(newdata[, c('gene_exp', 'LCI', 'UCI')])
      gene_explim <- range(gene_expcover[is.finite(gene_expcover) &
                                           gene_expcover != 0])
    }
  }
  df_long$f <- as.numeric(factor(df_long$f))
  maxf <- max(as.numeric(df_long$f), na.rm=T)
  df_long$time2 <- df_long$time + (df_long$f -1) * 2
  if (length(pch)==1) pch=rep(pch, maxf)
  if (length(bg)==1) bg=rep(bg, maxf)

  # Generate base plots
  if(graphics == "base"){
    # plot all samples
    if (plottype=="1" | class(fit)=='try-error') {
      plot(df_long$time2, df_long$gene_exp +.1,
           type='p', bty='l', las=1,
           xaxt='n', cex.axis=1.5, cex.lab=1.8,
           pch=pch[df_long$f], bg=bg[df_long$f],
           cex=1.8, xlab=NA, ylab=NA,
           log=switch(datatype, "TPM"="y", "VST"=""),
           #...,
           panel.first={
             for (i in paired) {lines(df_long$time2[df_long$id==i],
                                      df_long$gene_exp[df_long$id==i]+.1)}
           })
      #plot means
    } else {
      plot(newdata$time, newdata$gene_exp, type='p', bty='l', las=1,
           xaxt='n', cex.axis=1.5, cex.lab=1.8,
           pch=pch[newdata$f], bg=bg[newdata$f], cex=1.8, xlab=NA,
           ylim=gene_explim,
           log=switch(datatype, "TPM"="y", "VST"=""),
           panel.first={
             for (i in 1:maxf) {
               lines(newdata$time[newdata$f==i],
                     newdata$gene_exp[newdata$f==i],
                     col=bg[i], lwd=2)
             }
             arrows(x0=newdata$time, y0=newdata$LCI, y1=newdata$UCI,
                    angle=90, code=3, length=0.08)
           })
    }
    axis(1, 1:(2*maxf), labels=rep(time_labs, maxf), cex.axis=1.6)
    if (length(labs)==3) labs <- gsub('\\..*', '', labs)
    axis(1, 1:maxf*2-0.5, labels=labs, line=1.5, cex.axis=1.4, tick=F)
    if (xlab!="") mtext(xlab, 1, line=5, cex=1.2)
    if (ylab!="") mtext(ylab, 2, line=3, cex=1.4)
    if(title!="") mtext(title, side=3, adj=0, padj=-3, cex=1.5)
    if (suffix!="") suffix <- paste0(" [", suffix, "]")
    if (datatype!='VST' & is.na(dispersion)) {
      mtext("Dispersion NA", side=3, adj=0.04)
    } else if (class(fit)=='try-error') {
      mtext("Model error", side=3, adj=0.04)
    } else {
      mtext(bquote(paste("P"["time"]*"=", .(pval[1]), ", P"[.(plab)]*"=",
                         .(pval[2]), ", P"["time:"*.(plab)]*"=", .(pval[3]),
                         .(suffix))), side=3, adj=0.04)
    }

  } else{ # Generate ggplot plots
    if (plottype=="1" | class(fit)=='try-error') {
      p <- ggplot(df_long, aes(x=factor(time2), y=gene_exp + 0.1, group=id,
                               shape=factor(f), fill=factor(f),
                               color=factor(f))) +
        geom_point(size=3) +
        geom_line(color="black") +
        theme_classic() +
        scale_fill_manual(values=bg) +
        scale_color_manual(values=bg) +
        scale_shape_manual(values=pch) +
        labs(x=xlab, y=ylab, title=title) +
        theme(legend.position = "none") +
        scale_x_discrete(labels=setNames(rep(time_labs, 2), 1:4)) +
        geom_text(data=data.frame(label = labs, x = c(1.5, 3.5),
                                  y = min(df_long$gene_exp + 0.1)),
                  mapping=aes(label=label, x=x, y=y), hjust = 0.5,
                  vjust=5, nudge_x=0,  inherit.aes=F)  +
        coord_cartesian(clip = 'off')

      if(addBox) {
        p <- p + geom_boxplot(mapping=aes(x=time2, y=gene_exp + 0.1, group=time2),
                              color="black",  fill=NA)
      }


    } else{

      p <- ggplot(newdata, aes(x=time, y=gene_exp,
                               shape=factor(f), fill=factor(f),
                               color=factor(f))) +
        geom_point(size=3) +
        geom_line(color="black") +
        theme_classic() +
        scale_fill_manual(values=bg) + scale_color_manual(values=bg) +
        scale_shape_manual(values=pch) +
        labs(x=xlab, y=ylab, title=title) + theme(legend.position = "none") +
        scale_x_continuous(labels=setNames(rep(time_labs, 2), 1:4)) +
        geom_text(data=data.frame(label = labs, x = c(1.5, 3.5),
                                  y = min(df_long$gene_exp + 0.1)),
                  mapping=aes(label=label, x=x, y=y), hjust = 0.5,
                  vjust=5, nudge_x=0,  inherit.aes=F)  +
        coord_cartesian(clip = 'off') #+ lims(y=gene_explim)
    }

    if(datatype == "TPM") p <- p + scale_y_continuous(trans='log10')



    if (suffix!="") suffix <- paste0(" [", suffix, "]")
    if (datatype!='VST' & is.na(dispersion)) {
      p = p + labs(title= "Dispersion NA")
    } else if (class(fit)=='try-error') {
      p = p + labs(title= "Model error")
      mtext("Model error", side=3, adj=0.04)
    } else {
      p = p + labs(subtitle= bquote(
        paste("P"["time"]*"=", .(pval[1]),
              ", P"[.(plab)]*"=", .(pval[2]),
              ", P"["time:"*.(plab)]*"=",
              .(pval[3]), .(suffix))))
    }

    p

  }
}
