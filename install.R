#!/usr/bin/Rscript

# install.packages("devtools", repos = "http://cran.us.r-project.org")
# devtools::install()

install.packages(c('tidyr', 'dplyr', 'R.utils', 'rjson'), repos = 'https://cloud.r-project.org/')

# https://support.bioconductor.org/p/9150992/#9150993
BiocManager::install(version = "3.17")
BiocManager::install("EnrichmentBrowser")
