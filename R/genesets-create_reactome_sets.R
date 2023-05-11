#' create_reactome_set
#'
#' @param species string specifying the species to create the gene sets for
#'
#' @import dplyr
#' @import readr
#' @import glue
#' @export create_reactome_sets
create_reactome_sets <- function(species){
  reactome_url <- "https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt"

  species <- species_short_to_long(species)
  info(EnrichmentLogger(), glue("Download Reactome: {reactome_url} and subset species: {species}"))
  download.file(url = reactome_url, destfile = "UniProt2Reactome_All_Levels.txt")

  reactome_file <- read_delim("UniProt2Reactome_All_Levels.txt", col_names = FALSE)
  reactome_file <- reactome_file[reactome_file$X6 %in% species,]
  n_sets_species <- length(unique(reactome_file$X2))

  if(nrow(reactome_file) == 0){
    stop(glue("Species {species} not found in Reactome file."))
  }

  info(EnrichmentLogger(), glue("N: {n_sets_species} sets found for species: {species}"))

  colnames(reactome_file) <- c("items", "annotation_id","linkout", "annotation_name", "code", "organism")
  size_sets <- reactome_file |> dplyr::count(annotation_id, name = "size")

  reactome_file <- reactome_file |> dplyr::left_join(size_sets)
  reactome_file$date_download <- Sys.Date()

  gene_sets_table <- reactome_file |> dplyr::select(annotation_id, annotation_name, size, linkout, date_download) |> unique()
  genes_table <- reactome_file |> dplyr::select(annotation_id, items)


  if(nrow(gene_sets_table) != n_sets_species){
    stop(glue("Some sets have been lost in merging data."))
  }

  versions <- data.frame(annotation_name = "Reactome", download_link = reactome_url, version = Sys.Date())

  gene_sets_table$annotation_description <- NA
  reactomeSets <- new("AnnotationSetsReactome",
                       genesets_table = gene_sets_table,
                       genesets_proteins = genes_table,
                       versions = versions)

  return(reactomeSets)

}
