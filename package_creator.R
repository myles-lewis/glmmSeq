# Script to create the R package
# path_package = "/media/gcpeac/Katriona/seq_mixed_models/glmmSeq" # workstation
path_package = "~/Documents/Analyses/PhD_Projects/glmmSeq/" # laptop

library(devtools)
library(roxygen2)

## On inital set up

# remove any lingering objects from the environment
rm(fcPlot, aPlot, glmmQvals, glmmResults, glmmSeq, GlmmSeq, character_or_list, 
   df_or_matrix, maPlotly, modelPlot, labs, glmmGene, geneName, glmmResult, maPlot)
rm(list = c("predict"))
rm(list = c(".__C__GlmmSeq", "pairedPlot"))

# Set up documentation
devtools::document(path_package)
devtools::document(path_package)

# Build the package
devtools::build(path_package, vignettes=FALSE)

# Install
devtools::install(path_package)
devtools::build_vignettes(path_package)

library(glmmSeq)

# Check some of the documentation
?glmmSeq
?glmmQvals
?modelPlot
?glmmGene



# Test how "good" the package is
tests = goodpractice::all_checks()[! grepl("covr", goodpractice::all_checks())]
goodpractice::gp(checks=tests)
devtools::spell_check()
#spelling::update_wordlist(pkg = path_package)



# pkgdown
usethis::use_pkgdown()
pkgdown::build_site()


