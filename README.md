[![Lifecycle: Maturing](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-mediumpurple.svg)](https://choosealicense.com/licenses/mit/)
[![CRAN status](https://www.r-pkg.org/badges/version/glmmSeq)](https://cran.r-project.org/package=glmmSeq)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FKatrionaGoldmann%2FglmmSeq&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
[![GitHub issues](https://img.shields.io/github/issues/KatrionaGoldmann/glmmSeq.svg)](https://GitHub.com/KatrionaGoldmann/glmmSeq/issues/)
[![GitHub
tag](https://img.shields.io/github/tag/KatrionaGoldmann/glmmSeq.svg)](https://GitHub.com/KatrionaGoldmann/glmmSeq/tags/)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/glmmSeq?color=orange)](https://CRAN.R-project.org/package=glmmSeq)
[![Travis](https://img.shields.io/travis/KatrionaGoldmann/glmmSeq.svg)](https://github.com/KatrionaGoldmann/glmmSeq)

# glmmSeq 

<img src="https://katrionagoldmann.github.io/glmmSeq/logo.png" align="right" alt="" width="200" hspace="20" />



This R package is designed to model gene expression with a general linear mixed model (glmm). This allows us to include random effects as well as fixed effects. For the purpose of the package we use the `glmer` function from the [`lme4`](https://CRAN.R-project.org/package=lme4)
package which fits a glmm.

This package focuses in particular on changes in genes expression between different response or treatment groups over time. 


# Loading the package

### From CRAN

```
install.packages("glmmSeq")
```

### From Github

```
devtools::install_github("KatrionaGoldmann/glmmSeq")
```

### Locally

You can also download the source directory and load the functions individually:

```
functions = list.files("./R", full.names = TRUE)
invisible(lapply(functions, source))
```

But you will need to load in the additional libraries then:

```
# Install CRAN packages
invisible(lapply(c("MASS", "car", "ggplot2", "ggpubr", "lme4", 
                     "methods", "parallel", "plotly", "stats", 
                     "gghalves"),
                 function(p){
                   if(! p %in% rownames(installed.packages())) {
                     install.packages(p)
                   }
                   library(p, character.only=TRUE)
                 }))

# Install BioConductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
invisible(lapply(c("qvalue"), function(p){
  if(! p %in% rownames(installed.packages())) BiocManager::install(p)
  library(p, character.only=TRUE)
}))
```

# Example script

For examples see the [vignette](https://katrionagoldmann.github.io/glmmSeq/articles/glmmSeq.html). 

# Reference

glmmSeq was developed by the bioinformatics team at the [Experimental Medicine & Rheumatology department](https://www.qmul.ac.uk/whri/emr/) and [Centre for Translational Bioinformatics](https://www.qmul.ac.uk/c4tb/) at Queen Mary University London.

If you use this package please cite as:

```
citation("glmmSeq")

## To cite package ‘glmmSeq’ in publications use:
##
##  Myles Lewis, Katriona Goldmann, Elisabetta Sciacca, Cankut Cubuk and Anna Surace (2021). 
##  glmmSeq: General Linear Mixed Models for Gene-level Differential Expression. 
##  R package version 0.0.1. https://github.com/KatrionaGoldmann/glmmSeq
##
## A BibTeX entry for LaTeX users is
##
##  @Manual{,
##    title = {glmmSeq: General Linear Mixed Models for Gene-level Differential Expression},
##    author = {Myles Lewis and Katriona Goldmann and Elisabetta Sciacca and Cankut Cubuk and Anna Surace},
##    year = {2021},
##    note = {R package version 0.0.1},
##    url = {https://github.com/KatrionaGoldmann/glmmSeq},
##  }
```

