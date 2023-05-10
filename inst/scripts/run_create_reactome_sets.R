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

	make_option(c("-s","--reactome_species"), type = "character", default = NULL,
		help = "Character specifying the species to create the sets for.")

);

# Parse arguments
opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);

print(glue("Species: {opt$reactome_species}"))

date <- Sys.Date()
reactome_sets <- EnrichmentAnalysisStepR::create_reactome_sets(opt$reactome_species)

taxo <- species_to_taxonomy(opt$reactome_species)
outdir <- glue("inst/extdata/{taxo}/")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save(reactome_sets, file = glue("{outdir}/Reactome_{date}.rda"))

EnrichmentAnalysisStepR::write_gmt(reactome_sets$genesets_proteins, outfolder = outdir, db_name = "Reactome")
write_set_info_reactome(reactome_sets = reactome_sets, outfolder = outdir)


