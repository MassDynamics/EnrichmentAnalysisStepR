#' Taxonomy to long species name
#' @param taxonomy taxonomy code
#' @export
taxonomy_to_long_species <- function(taxonomy){
  species <- switch(
    taxonomy,
    "9606" = "Homo sapiens",
    "10090" = "Mus musculus",
    "559292" = "Yeast",
    "10029" = "Chinese hamster"
  )
  return(species)
}

#SPECIES = ['Human', 'Mouse', 'Yeast', 'Chinese hamster', SPECIES_OTHER].freeze
#TAXON_IDS = %w[9606 10090 559292 10029 9606].freeze


#' Taxonomy to short species name
#' @param taxonomy taxonomy code
#' @export
taxonomy_to_short_species <- function(taxonomy){
  species <- switch(
    taxonomy,
    "9606" = "Human",
    "10090" = "Mouse",
    "559292" = "Yeast",
    "10029" = "Chinese hamster"
  )
  return(species)
}


#' Species to taxonomy
#' @param species taxonomy code
#' @export
species_to_taxonomy <- function(species){
  taxonomy <- switch(
    species,
    "Human" = "9606",
    "Homo sapiens" = "9606",
    "human" = "9606",
    "Mouse" = "10090",
    "mouse" = "10090",
    "Yeast" = "559292",
    "yeast" = "559292",
    "Chinese hamster" = "10029"
  )
  return(taxonomy)
}

#' Species short to long names
#' @param taxonomy taxonomy code
#' @export
species_short_to_long <- function(taxonomy){
  species <- switch(
    taxonomy,
    "Human" = "Homo sapiens",
    "human" = "Homo sapiens",
    "Mouse" = "Mus musculus",
    "mouse" = "Mus musculus",
    "Yeast" = "Saccharomyces cerevisiae",
    "yeast" = "Saccharomyces cerevisiae",
    "Chinese hamster" = NULL
  )
  return(species)
}

#' Species to GAF
#' @param species species
#' @export
species_to_gaf <- function(species){
  gaf_id <- switch(
    species,
    "Human" = "goa_human",
    "human" = "goa_human",
    "Mouse" = "mgi",
    "mouse" = "mgi", # Mus musculus
    "yeast" = "sgd", # http://current.geneontology.org/annotations/sgd.gaf.gz
    "Chinese hamster" = NULL
  )
  return(gaf_id)
}


#' Species to GAF URL
#' @param species species
#' @export
species_to_url <- function(species){
  url <- switch(
    species,
    "Human" = "http://geneontology.org/gene-associations/goa_human.gaf.gz",
    "human" = "http://geneontology.org/gene-associations/goa_human.gaf.gz",
    "Mouse" = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz",
    "mouse" = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/MOUSE/goa_mouse.gaf.gz", # Mus musculus
    "yeast" = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/YEAST/goa_yeast.gaf.gz", # http://current.geneontology.org/annotations/sgd.gaf.gz
    "Chinese hamster" = NULL
  )
  return(url)
}

# https://chomine.boku.ac.at/chomine/dataCategories.do
