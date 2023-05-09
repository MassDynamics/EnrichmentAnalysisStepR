#!/bin/bash

Rscript -e "devtools::build(binary=FALSE); install.packages('/Users/annaquaglieri/Projects/MD/MD-services/EnrichmentAnalysisStepR_0.0.26.tar.gz', repos = NULL, type='source')"

echo "MsigDB Human v7.5"
make msigdb MSIG_VERSION="7.5"

echo "Go sets"
make go SPECIES="human"
make go SPECIES="mouse"
make go SPECIES="yeast"

echo "Reactome sets"
make reactome SPECIES="human"
make reactome SPECIES="mouse"
make reactome SPECIES="yeast"


