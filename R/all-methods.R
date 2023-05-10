# Create a generic getConditionComparisons function
setGeneric("getConditionComparisons", function(params) {
  standardGeneric("getConditionComparisons")
})

setGeneric("getGeneSets", function(params) {
  standardGeneric("getGeneSets")
})

setGeneric("getSpecies", function(params) {
  standardGeneric("getSpecies")
})


setGeneric("subsetCategory", function(.Object, category) {
  standardGeneric("subsetCategory")
})

setGeneric("getIntensities", function(dataset) {
  standardGeneric("getIntensities")})

#' getConditionComparisons method
#' @import dplyr
#' @export
setMethod("getConditionComparisons", signature(params = "EnrichmentDatasetParameters"), function(params){
  return(dplyr::bind_rows(params@fields$condition_comparisons))
})

#' getGeneSets method
#' @export
setMethod("getGeneSets", signature(params = "EnrichmentDatasetParameters"), function(params){
  sets <- params@fields$species_gene_set_selection$gene_sets
  sets <- do.call(c, sets)
  return(sets)
})


#' getSpecies method
#' @export
setMethod("getSpecies", signature(params = "EnrichmentDatasetParameters"), function(params){
  return(params@fields$species_gene_set_selection$species)
})


#' subsetCategory in gene sets method
#' @export
setMethod("subsetCategory", "AnnotationSetsMsigDB", function(.Object, category){
  newObject <- .Object
  newObject@genesets_table <- newObject@genesets_table[newObject@genesets_table$set_category %in% category,]
  subset_sets_ids <- unique(newObject@genesets_table$annotation_id)
  newObject@genesets_proteins <- newObject@genesets_proteins[newObject@genesets_proteins$annotation_id %in% subset_sets_ids,]
  return(newObject)
})


#' Get Intensities data from IntensitiesTable class
#' @export
setMethod("getIntensities", signature(dataset = "IntensitiesTable"), function(dataset){
  dataset@data
})

