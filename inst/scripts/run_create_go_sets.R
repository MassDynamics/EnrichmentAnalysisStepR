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

	make_option(c("-s","--species"), type = "character", default = NULL,
		help = "Url of GO gab file to download, depending on the species.
		See http://current.geneontology.org/products/pages/downloads.html."),

	make_option(c("-t","--save_txt_files"), type = "logical", default = FALSE,
	            help = "Save GMT and info.csv file like in the old Mass Dynamics annotation format.")

);


# http://geneontology.org/gene-associations/goa_human.gaf.gz

# Parse arguments
opt_parser = OptionParser(option_list=option_list, add_help_option = TRUE);
opt = parse_args(opt_parser);

date <- Sys.Date()
species <- opt$species
goSets <- EnrichmentAnalysisStepR::create_go_sets(species = species)

taxo <- species_to_taxonomy(species)
outdir <- glue("inst/extdata/{date}/{taxo}/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
save(goSets, file = glue("{outdir}/GO.rda"))

categories <- names(table(goSets@genesets_table$annotation_category))
for(category in categories){
  subsetSet <- subsetCategory(goSets, category = category)
  category <- EnrichmentAnalysisStepR:::mapGOCategoryNames(category)
  save(subsetSet, file = glue("{outdir}/GO_{category}.rda"))
}

# Save old format files
if(opt$save_txt_files){
  write_gmt_go(go_sets = goSets, outfolder = outdir)
  write_set_info_go(go_sets = goSets, outfolder = outdir)
}




