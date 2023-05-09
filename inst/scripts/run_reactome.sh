#!/bin/bash

species=$1

echo $species

# Rscript -e "devtools::build(binary=FALSE); install.packages('/Users/annaquaglieri/Projects/MD/MD-services/EnrichmentAnalysisStepR_0.0.26.tar.gz', repos = NULL, type='source')"

Rscript ./inst/scripts/run_create_reactome_sets.R --reactome_species ${species}
