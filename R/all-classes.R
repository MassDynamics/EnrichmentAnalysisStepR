# http://adv-r.had.co.nz/S4.html

#' EnrichmentDatasetParameters class
#' @export
setClass("EnrichmentDatasetParameters", slots = list(dataset_id = "character", version = "character", fields = "list"))


setMethod("initialize", "EnrichmentDatasetParameters", function(.Object,
                                                                dataset_id,
                                                                version,
                                                                fields) {
  # Define the required column names
  field_names <- c("experiment_design", "species_gene_set_selection", "condition_comparisons")

  # Check if the input data frame has the required column names
  if (!all(field_names %in% names(fields))) {
    stop("Fields in dataset parameters should have the following fields: ", paste(field_names, collapse = ", "))
  }

  .Object@dataset_id <- dataset_id
  .Object@version <- version
  .Object@fields <- fields

  return(.Object)
})


#' AnnotationSets class
#' @export
setClass(
  "AnnotationSets",
  slots = list(genesets_table = "data.frame",
               genesets_proteins = "data.frame",
               versions = "data.frame")
)

#' AnnotationSets class
#' @export
setClass(
  "AnnotationSetsMsigDB",
  slots = list(genes_protein_mapping = "data.frame"), contains = "AnnotationSets"
)

setMethod("initialize", "AnnotationSetsMsigDB", function(.Object, genesets_table,
                                                         genesets_proteins,
                                                         versions,
                                                         genes_protein_mapping) {
  # Define the required column names
  required_genesets_table <- c("annotation_name", "annotation_id", "organism", "set_category", "set_subcategory", "set_exact_source",
                        "linkout", "annotation_description", "annotation_description_full", "size_initial_set", "size", "prop_unmapped_items" )

  # Check if the input data frame has the required column names
  if (!all(required_genesets_table %in% colnames(genesets_table))) {
    stop("genesets_table must have the following column names: ", paste(required_genesets_table, collapse = ", "))
  }

  # Define the required column names
  required_genesets_proteins <- c("annotation_id", "items")
  if (!all(required_genesets_proteins %in% colnames(genesets_proteins))) {
    stop("genesets_proteins must have the following column names: ", paste(required_genesets_proteins, collapse = ", "))
  }

  # Assign the data frame to the data slot
  .Object@genesets_table <- genesets_table
  .Object@genesets_proteins <- genesets_proteins
  .Object@versions <- versions
  .Object@genes_protein_mapping <- genes_protein_mapping

  return(.Object)
})


#######################
##### IntensitiesTable
#' @export
setClass("IntensitiesTable", slots = c(data = "data.frame"))




