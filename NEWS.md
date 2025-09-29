News
=====

# glmmSeq 0.5.7
###### 29/09/2025
* Fixing broken links in documetation.

# glmmSeq v0.5.6
###### 17/01/2025
* Now returns variance of the random effect term as an additional column in the 
`@predict` slot.

# glmmSeq v0.5.5
###### 29/09/2022
* Improve error handling if all genes give errors

# glmmSeq v0.5.4
###### 25/09/2022
* Add option for likelihood ratio test (LRT) in `glmmSeq()` and `lmmSeq()`
* Enable changing of model formula, family and control parameters when refitting
model for a single gene in `glmmRefit()`

# glmmSeq v0.5.2
###### 19/09/2022
* Fix passing of `...` in parallelisation on windows in `glmmSeq()` & `lmmSeq()`

# glmmSeq v0.5.1
###### 02/09/2022
* Add option to `glmmSeq()` to use `glmmTMB` package for fitting negative
binomial GLMM (or other GLM family) models
* Add `plab` argument to `modelPlots` to customise p-value labels
* Add mean expression column `meanExp` to results in `@stats$res` slot

# glmmSeq v0.4.0
###### 10/08/2022
* Further speed enhancements to `lmmSeq()` using `lme4::modular` code. Speed
increase of around 25%.

# glmmSeq v0.3.0
###### 05/08/2022
* Significant update including output of standard error on fitted coefficients
for both `glmmSeq` and `lmmSeq`
* Add `summary` function for glmmSeq and lmmSeq to display results for an
individual gene
* Reorganised `@stats` slot in output objects to include more information
including DF
* Added option of using Saiterthwaite's DF method with ANOVA type III tables as
option for `lmmSeq` using the `lmerTest` package
* Added `lmmRefit` function to fit an identical (g)lmer model. This can then be
passed to the `emmeans` package for visualisation of more complex models.
* Separated `modelPlots` (base graphics) and `ggmodelPlots` (ggplot2)
* Streamlined `modelPlots` to allow for simplest case `gene ~ Time + (1 | ID)`

# glmmSeq v0.2.2
###### 16/07/2022
* Fast version of `lmmSeq()`. Improves speed of calculation of type 2 Wald test.
Overall speed increase of 25-50%.
* Faster version of `glmmSeq()`
* Automatically detects `id` column name from the RE term in the formula 

# glmmSeq v0.2.1
###### 11/07/2022
* Add `...` option to `fcPlot` which is passed to `plotly()` or `ggplot()`
* set `annotationPosition=FALSE` in `fcPlot` so arrows/connectors are not moved

# glmmSeq v0.2.0
###### 09/07/2022
* Add lmmSeq function for gaussian linear mixed models

# glmmSeq v0.1.2
###### 08/07/2022
* Fix CMD check errors

# glmmSeq v0.1.1
###### 18/04/2021

##### Bug Fixes
* Fix colours in pairedPlots and modelPlots

# glmmSeq v0.1.0
###### 16/03/2021

##### Features
* Add progress bars to glmmSeq functions
* Add option `returnList` to return glmmSeq output as list (to make error
catching easier)

##### Bug Fixes
* Fix missing plots in vignette

# glmmSeq v0.0.1
###### 05/03/2021

* This is the initial build of `glmmSeq`
