#!/usr/bin/env Rscript

install.packages("optparse", repos='http://cran.us.r-project.org')
library(optparse)
library(PrepareGeneSets)
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
msigdb_version <- gsub(".xml", "", basename(msigdb_file))

msig_sets <- PrepareGeneSets::generate_msig_sets(msigdb_file, msigdb_version, opt$uniprot_table_path)

version <- opt$mdigdb_version
species <- "human"
taxo <- species_to_taxonomy(species)
outdir <- glue("inst/extdata/{taxo}/")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save(msig_sets, file = glue("{outdir}/{msigdb_version}.rda"))

