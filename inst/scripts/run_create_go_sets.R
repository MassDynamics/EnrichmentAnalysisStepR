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

	make_option(c("-s","--species"), type = "character", default = NULL, 
		help = "Url of GO gab file to download, depending on the species. 
		See http://current.geneontology.org/products/pages/downloads.html.")
	
); 


# http://geneontology.org/gene-associations/goa_human.gaf.gz

# Parse arguments
opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);

date <- Sys.Date()
species <- opt$species
go_sets <- PrepareGeneSets::create_go_sets(species = species)

taxo <- species_to_taxonomy(species)
outdir <- glue("inst/extdata/{taxo}/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save(go_sets, file = glue("{outdir}/GO_{species}_{date}.rda"))

write_gmt_go(go_sets = go_sets, outfolder = outdir)
write_set_info_go(go_sets = go_sets, outfolder = outdir)
