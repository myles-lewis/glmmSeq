# Script to create the R package
 # Set wd to pacakge directory

# restart session
.rs.restartR()

library(devtools)
library(roxygen2)

## On inital set up

# remove any lingering objects from the environment
rm(list=ls(all=T))

# Set up documentation
devtools::document()
devtools::document()

# Build the package
devtools::build(vignettes=FALSE)

# Install
devtools::install()
devtools::build_vignettes()

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
spelling::update_wordlist()



# pkgdown
pkgdown::build_site()


