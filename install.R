#!/usr/bin/Rscript

# install.packages("devtools", repos = "http://cran.us.r-project.org")
# devtools::install()

install.packages(c('tidyr', 'dplyr', 'R.utils', 'rjson'), repos = 'https://cloud.r-project.org/')
BiocManager::install("EnrichmentBrowser")
