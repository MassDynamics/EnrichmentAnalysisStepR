#!/usr/bin/env Rscript

install.packages("optparse", repos='http://cran.us.r-project.org')
library(optparse)
library(EnrichmentAnalysisStepR)
library(glue)

####################
# Read in arguments
####################

print("Reading arguments in...")

# Useful when a matched normal is not available. In this way the variants called will be annotated with the PON.
# Define options for command line
option_list = list(

	make_option(c("-p","--xml_filepath"), type = "character", default = NULL,
		help = "Path to the XML MsigDB file. From: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"),

	make_option(c("-v","--msigdb_version"), type = "character", default = NULL,
	            help = "msigdb version"),

  make_option(c("-u", "--uniprot_table_path"), type = "character", default = NULL,
              help = "Full path to uniprot table.")

);

# Parse arguments
opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);

msigdb_file <- opt$xml_filepath
msigdb_filename <- gsub(".xml", "", basename(msigdb_file))

msigSets <- EnrichmentAnalysisStepR::generate_msig_sets(msigdb_file, msigdb_filename, opt$uniprot_table_path)

version <- opt$msigdb_version
species <- switch(version,
                  "7.5" = "human",
                  "2023.1.Mm" = "mouse")
taxo <- species_to_taxonomy(species)

date <- Sys.Date()
outdir <- glue("inst/extdata/{date}/{taxo}/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save(msigSets, file = glue("{outdir}/{msigdb_filename}.rda"))

categories <- names(table(msigSets@genesets_table$annotation_category))
for(category in categories){
  subsetSet <- subsetCategory(msigSets, category = category)
  save(subsetSet, file = glue("{outdir}/{msigdb_filename}_{category}.rda"))
}

