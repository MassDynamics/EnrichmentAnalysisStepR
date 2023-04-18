datasets <- function() {}

#' Gets the condition comparisons provided by the dataset parameters
#' @param datasetParameters a jsonlite representation of the dataset parameters
#' @return a data frame representing the condition comparisons
#' @export datasets.getConditionComparisons
datasets.getConditionComparisons <- function(datasetParameters) {
  return(datasetParameters$fields$condition_comparisons)
}

#' Gets the gene sets information provided by the dataset parameters
#' @param datasetParameters a jsonlite representation of the dataset parameters
#' @return a list representing the gene sets information
#' @export datasets.getGeneSetsInfo
datasets.getGeneSetsInfo <- function(datasetParameters) {
  return(
    list(
      species = datasetParameters$fields$species_gene_set_selection$species,
      gene_set_names = datasetParameters$fields$species_gene_set_selection$gene_sets
    )
  )
}