#!/bin/bash

set -e

build=$1
msig_version=$2

if [ "$build" == "true" ]; then
  Rscript -e "devtools::build(binary=FALSE); install.packages('/Users/annaquaglieri/Projects/MD/MD-services/EnrichmentAnalysisStepR_0.0.26.tar.gz', repos = NULL, type='source')"
fi

echo "#########"
echo "Create GO"
echo "#########"

Rscript ./inst/scripts/run_create_go_sets.R --species "human"
Rscript ./inst/scripts/run_create_go_sets.R --species "mouse"
Rscript ./inst/scripts/run_create_go_sets.R --species "yeast"

echo "#########"
echo "Create Reactome"
echo "#########"

Rscript ./inst/scripts/run_create_reactome_sets.R --reactome_species "human"
Rscript ./inst/scripts/run_create_reactome_sets.R --reactome_species "mouse"
Rscript ./inst/scripts/run_create_reactome_sets.R --reactome_species "yeast"


echo "#########"
echo "Create MsigdbDB"
echo "#########"

echo "Download Uniprot"
aws s3 cp s3://md-external-knowledge/uniprot-versions/uniprot-2Creviewed_2C-2023.01.19-04.19.19.92.tsv .
uniprot_table_path=uniprot-2Creviewed_2C-2023.01.19-04.19.19.92.tsv

msigdb_xml_path=msigdb_v${msig_version}.xml
if [ ! -e "$msigdb_xml_path" ]; then
  echo "Download MsigDB"
  wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/${msig_version}/msigdb_v${msig_version}.xml -O msigdb_v${msig_version}.xml
fi
echo $msigdb_xml_path

Rscript ./inst/scripts/run_create_msigdb_sets.R \
--xml_filepath ${msigdb_xml_path} \
--msigdb_version ${msig_version} \
--uniprot_table_path ${uniprot_table_path}
