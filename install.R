#!/usr/bin/Rscript

# install.packages("devtools", repos = "http://cran.us.r-project.org")
# devtools::install()

install.packages(c('tidyr', 'dplyr', 'R.utils', 'rjson'), repos = 'https://cloud.r-project.org/')

# We need to upgrade R to 4.3
# https://support.bioconductor.org/p/9150992/#9150993
# BiocManager::install(version = "3.17")
# BiocManager::install("EnrichmentBrowser")

# See older version install https://www.bioconductor.org/install/
source("https://bioconductor.org/biocLite.R")
BiocInstaller::biocLite(c("EnrichmentBrowser"))
