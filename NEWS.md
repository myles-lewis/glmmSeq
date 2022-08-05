News
=====

# glmmSeq v0.3.0
###### 05/8/2022
* Significant update including output of standard error on fitted coefficients
for both `glmmSeq` and `lmmSeq`
* Add `summary` function for glmmSeq and lmmSeq to display results for an
individual gene
* Reorganised `@stats` slot in output objects to include more information
including DF
* Added option of using Saiterthwaite's DF method with anova type III tables as
option for `lmmSeq` using the `lmerTest` package
* Added `lmmRefit` function to fit an identical (g)lmer model. This can then be
passed to the `emmeans` package for visualisation of more complex models.
* Separated `modelPlots` (base graphics) and `ggmodelPlots` (ggplot2)
* Streamlined `modelPlots` to allow for simplest case gene ~ Time + (1 | ID)

# glmmSeq v0.2.2
###### 16/7/2022
* Fast version of `lmmSeq()`. Improves speed of calculation of type 2 Wald test.
Overall speed increase of 25-50%.
* Faster version of `glmmSeq()`
* Automatically detects `id` column name from the RE term in the formula 

# glmmSeq v0.2.1
###### 11/7/2022
* Add `...` option to fcPlot which is passed to `plotly()` or `ggplot()`
* set `annotationPosition=FALSE` in `fcPlot so arrows/connectors are not moved

# glmmSeq v0.2.0
###### 09/7/2022
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
