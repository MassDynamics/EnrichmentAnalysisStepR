#!/bin/bash

species=$1

echo $species

#Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org'); devtools::build(binary=FALSE); install.packages('/Users/annaquaglieri/Projects/MD/MD-services/PrepareGeneSets_0.1.3.tar.gz', repos = NULL, type='source')"

Rscript ./inst/scripts/run_create_reactome_sets.R --reactome_species ${species}