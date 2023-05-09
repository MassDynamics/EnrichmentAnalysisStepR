#!/bin/bash

# http://geneontology.org/gene-associations/goa_human.gaf.gz
species=$1

echo $species

# Rscript -e "devtools::build(binary=FALSE); install.packages('/Users/annaquaglieri/Projects/MD/MD-services/EnrichmentAnalysisStepR_0.0.26.tar.gz', repos = NULL, type='source')"

Rscript ./inst/scripts/run_create_go_sets.R --species ${species}
