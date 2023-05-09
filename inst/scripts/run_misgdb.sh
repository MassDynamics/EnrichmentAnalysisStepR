#!/bin/bash

msig_version=$1

aws s3 cp s3://md-external-knowledge/uniprot-versions/uniprot-2Creviewed_2C-2023.01.19-04.19.19.92.tsv . 
uniprot_table_path=uniprot-2Creviewed_2C-2023.01.19-04.19.19.92.tsv

echo "MsigdbDB download"
echo $msig_version
wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/${msig_version}/msigdb_v${msig_version}.xml -O msigdb_v${msig_version}.xml
msigdb_xml_path=msigdb_v${msig_version}.xml

#Rscript -e "install.packages('devtools', repos = 'http://cran.us.r-project.org'); devtools::build(binary=FALSE); install.packages('/Users/annaquaglieri/Projects/MD/MD-services/PrepareGeneSets_0.1.3.tar.gz', repos = NULL, type='source')"

Rscript ./inst/scripts/run_create_msigdb_sets.R --xml_filepath ${msigdb_xml_path} --msigdb_version ${msig_version} --uniprot_table_path ${uniprot_table_path}
