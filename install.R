#!/usr/bin/Rscript

install.packages("renv", repos = "http://cran.us.r-project.org")
renv::install("bioc::Biobase", "bioc::EnrichmentBrowser")

renv::install()
