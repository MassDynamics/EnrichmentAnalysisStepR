# production environment install script

cat("Preparing R enviroment")

list.of.packages <- c(
  "jsonlite",
  "tidyr",
  "dplyr",
  "stringr",
  "R.utils",
  "BiocManager"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")


BiocManager::install("EnrichmentBrowser")
