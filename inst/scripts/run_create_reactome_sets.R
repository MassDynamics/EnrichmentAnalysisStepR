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
		help = "Character specifying the species to create the sets for."),

	make_option(c("-t","--save_txt_files"), type = "logical", default = FALSE,
	            help = "Save GMT and info.csv file like in the old Mass Dynamics annotation format.")

);

# Parse arguments
opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);
date <- Sys.Date()
species <- opt$reactome_species

print(glue("Species: {opt$reactome_species}"))
reactomeSets <- EnrichmentAnalysisStepR::create_reactome_sets(species)

taxo <- species_to_taxonomy(species)
outdir <- glue("inst/extdata/{taxo}/{date}")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save(reactomeSets, file = glue("{outdir}/Reactome.rda"))

EnrichmentAnalysisStepR::write_gmt(reactomeSets@genesets_proteins, outfolder = outdir, db_name = "Reactome")
save(reactomeSets, file = glue("{outdir}/Reactome.rda"))

# Save old format files
if(opt$save_txt_files){
  write_gmt(genesets_proteins = reactomeSets@genesets_proteins, outfolder = outdir, db_name = "Reactome")
  write_set_info_reactome(reactome_sets = reactomeSets, outfolder = outdir)
}


