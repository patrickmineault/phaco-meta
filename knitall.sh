#!/bin/sh
Rscript -e "rmarkdown::render('investigate-washout.Rmd', output_format='all')"
Rscript -e "rmarkdown::render('methods.Rmd', output_format='all')"
Rscript -e "rmarkdown::render('metaanalysis.Rmd', output_format='all')"

