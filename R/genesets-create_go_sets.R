#' create_go_set
#'
#' @param go_url url to gaf file
#'
#' @import dplyr
#' @import readr
#' @import ontologyIndex
#' @export create_go_sets
create_go_sets <- function(species = "human"){
  # http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
  # obo file for annotation descriptions http://geneontology.org/docs/download-ontology/ - use basic obo

  species <- tolower(species)
  go_url <- species_to_url(species)
  log4r::info(EnrichmentLogger(), glue("Download GO annotation: {go_url}"))
  dest_name <- basename(go_url)
  download.file(url = go_url, destfile = dest_name)

  log4r::info(EnrichmentLogger(), "Get GO sets information from OBO file")
  sets_information <- get_sets_information()

  log4r::info(EnrichmentLogger(), glue("Read GAF file for species: {species}"))
  x <- readLines(dest_name, n = 50)
  #start_line_pattern <- species_to_obo(species)
  start_line <- min(grep(paste0('^UniProt'), x))

  go <- read_delim(dest_name, delim = "\t",
                         escape_double = FALSE, col_names = FALSE,
                         trim_ws = TRUE, skip = start_line-1)

  go <- go[,c("X2","X5")]
  colnames(go) <- c("items", "annotation_id")
  go <- go %>% dplyr::left_join(sets_information)

  size_sets <- go |> dplyr::count(annotation_id, name = "size")
  go <- go |> dplyr::left_join(size_sets)
  go$date_download <- Sys.Date()
  go$linkout <- paste0("http://amigo.geneontology.org/amigo/term/",go$annotation_id)

  genesets_table <- go |> dplyr::select(annotation_id, annotation_name, set_subcategory, annotation_description, size, linkout, date_download) |> unique()
  genesets_table <- genesets_table %>% dplyr::rename(annotation_subcategory = set_subcategory)

  if(length(unique(go$annotation_id)) != nrow(genesets_table)){
    stop("Some GO annotations have been lost when creating the final data.")
  }

  genes_table <- go |> dplyr::select(annotation_id, items)

  versions <- data.frame(annotation_name = "GO", download_link = go_url, version = Sys.Date())

  list(genesets_table = genesets_table,
       genesets_proteins = genes_table,
       versions = versions)

}

#' Get sets information
#'
#' @param go_url url to gaf file
#'
#' @importFrom purrr reduce
#' @import ontologyIndex
#' @export get_sets_information
get_sets_information <- function(){
  obo_file_url <- "http://current.geneontology.org/ontology/go-basic.obo"
  log4r::info(EnrichmentLogger(), glue("OBO file from: {obo_file_url}"))
  dest_file <- "go-basic.obo"
  download.file(obo_file_url, destfile = dest_file)
  log4r::info(EnrichmentLogger(), glue("Check that downloaded OBO file exists: {file.exists(dest_file)}"))

  obo_file <- get_ontology("go-basic.obo", extract_tags = "everything")
  categories <- get_sets_category(obo_file)
  ann_names <- get_sets_names(obo_file)
  des <- get_sets_descriptions(obo_file)
  output <- purrr::reduce(list(categories,ann_names,des), dplyr::left_join, by = 'annotation_id')
  return(output)
}

#' Get GO sets category
#'
#' @param obo_file obo file imported with `get_ontology`
#' @export get_sets_category
get_sets_category <- function(obo_file){
  namespaces <- obo_file$namespace
  annotation_ids <- names(namespaces)
  category <- do.call(c, namespaces)
  return(data.frame(annotation_id = annotation_ids, set_subcategory = category))
}

#' Get GO sets descriptions
#'
#' @param obo_file obo file imported with `get_ontology`
#' @export get_sets_descriptions
get_sets_descriptions <- function(obo_file){
  descriptions <- obo_file$def
  descriptions <- gsub('\"', '', descriptions, fixed=TRUE)
  annotation_ids <- names(descriptions)
  return(data.frame(annotation_id = annotation_ids, annotation_description = descriptions))
}

#' Get GO sets names
#'
#' @param obo_file obo file imported with `get_ontology`
#' @export get_sets_names
get_sets_names <- function(obo_file){
  names <- obo_file$name
  annotation_ids <- names(names)
  return(data.frame(annotation_id = annotation_ids, annotation_name = names))
}




