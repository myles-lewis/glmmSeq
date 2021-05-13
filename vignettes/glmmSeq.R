## ----setup, include=FALSE---------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 8, fig.height = 6)
options(width=96)
library(kableExtra)

## ---- eval=FALSE------------------------------------------------------------------------------
#  install.packages("glmmSeq")

## ---- eval=FALSE------------------------------------------------------------------------------
#  devtools::install_github("KatrionaGoldmann/glmmSeq")

## ---- eval=FALSE------------------------------------------------------------------------------
#  functions = list.files("./R", full.names = TRUE)
#  invisible(lapply(functions, source))

## ---- eval=FALSE------------------------------------------------------------------------------
#  # Install CRAN packages
#  invisible(lapply(c("MASS", "car", "ggplot2", "ggpubr", "lme4", "methods",
#                     "parallel", "plotly", "stats", "gghalves"),
#                   function(p){
#                     if(! p %in% rownames(installed.packages())) {
#                       install.packages(p)
#                     }
#                     library(p, character.only=TRUE)
#                   }))
#  
#  # Install BioConductor packages
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  invisible(lapply(c("qvalue"), function(p){
#    if(! p %in% rownames(installed.packages())) BiocManager::install(p)
#    library(p, character.only=TRUE)
#  }))
#  

## ---- message=FALSE, warning=FALSE------------------------------------------------------------
library(glmmSeq)
set.seed(1234)

## ---------------------------------------------------------------------------------------------
data(PEAC_minimal_load)

## ---------------------------------------------------------------------------------------------
metadata$EULAR_binary  = NA
metadata$EULAR_binary[metadata$EULAR_6m %in%
                        c("Good responder", "Moderate responder" )] = "responder"
metadata$EULAR_binary[metadata$EULAR_6m %in% c("Non responder")] = "non_responder"
metadata = metadata[! is.na(metadata$EULAR_binary), ]

kable(head(metadata), row.names = F) %>% kable_styling()

## ---------------------------------------------------------------------------------------------
tpm = tpm[, metadata$SAMID]
kable(head(tpm)) %>% kable_styling() %>%
  scroll_box(width = "100%")

## ---------------------------------------------------------------------------------------------
disp <- apply(tpm, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
  })

head(disp)

## ---- message=FALSE---------------------------------------------------------------------------
disp  <- setNames(edgeR::estimateDisp(tpm)$tagwise.dispersion, rownames(tpm))

head(disp)

## ---- eval=FALSE------------------------------------------------------------------------------
#  dds <- DESeqDataSetFromTximport(txi = txi, colData = metadata, design = ~ 1)
#  dds <- DESeq(dds)
#  dispersions <- setNames(dispersions(dds), rownames(txi$counts))

## ---------------------------------------------------------------------------------------------
sizeFactors <- colSums(tpm)  
sizeFactors <- sizeFactors / mean(sizeFactors)  # normalise

head(sizeFactors)

## ---- eval=FALSE------------------------------------------------------------------------------
#  sizeFactors <- calcNormFactors(counts, method="TMM")

## ---- eval=FALSE------------------------------------------------------------------------------
#  sizeFactors <- estimateSizeFactorsForMatrix(counts)

## ---- warning=FALSE---------------------------------------------------------------------------
results <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                  id = "PATID",
                  countdata = tpm,
                  metadata = metadata,
                  dispersion = disp,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = 1)

## ---- warning=FALSE---------------------------------------------------------------------------
results2 <- glmmSeq(~ Timepoint * EULAR_binary + (1 | PATID),
                  id = "PATID",
                  countdata = tpm,
                  metadata = metadata,
                  dispersion = disp,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  cores = 1)

## ---------------------------------------------------------------------------------------------
names(attributes(results))

## ---------------------------------------------------------------------------------------------
kable(results@modelData) %>% kable_styling()

## ---------------------------------------------------------------------------------------------
stats = data.frame(results@stats)

kable(stats[order(stats$P_Timepoint.EULAR_6m), ]) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

## ---------------------------------------------------------------------------------------------
predict = data.frame(results@predict)
kable(predict) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "400px")

## ---------------------------------------------------------------------------------------------
results <- glmmQvals(results, pi0=1)

## ---- warning=FALSE---------------------------------------------------------------------------
MS4A1glmm <- glmmSeq(~ Timepoint * EULAR_6m + (1 | PATID),
                     id = "PATID",
                     countdata = tpm["MS4A1", ],
                     metadata = metadata,
                     dispersion = disp,
                     verbose=FALSE)

## ---- warning=FALSE---------------------------------------------------------------------------
MS4A1fit <- glmmGene(~ Timepoint * EULAR_6m + (1 | PATID),
                     gene = "MS4A1",
                     id = "PATID",
                     countdata = tpm,
                     metadata = metadata,
                     dispersion = disp['MS4A1'])

MS4A1fit

