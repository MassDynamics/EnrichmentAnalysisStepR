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

#' createEnrichmentDatasetParameters
#' @export
createEnrichmentDatasetParameters <- function(dataset_id, version, fields) {
  return(new("EnrichmentDatasetParameters", dataset_id = dataset_id, version = version, fields = fields))
}


########################
#' AnnotationSets class
########################
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
  required_genesets_table <- c(getSharedColNamesAnnotations(),
                               "organism", "annotation_category",
                               "annotation_subcategory", "annotation_exact_source",
                               "annotation_description_full", "size_initial_set","prop_unmapped_items")

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


#' AnnotationSetsGo class
#' @export
setClass(
  "AnnotationSetsGO", contains = "AnnotationSets"
)

setMethod("initialize", "AnnotationSetsGO", function(.Object, genesets_table,
                                                         genesets_proteins,
                                                         versions) {

  # Define the required column names
  required_genesets_table <- c(getSharedColNamesAnnotations(), "annotation_category", "date_download" )

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

  return(.Object)
})

#' AnnotationSetsGo class
#' @export
setClass(
  "AnnotationSetsReactome", contains = "AnnotationSets"
)

setMethod("initialize", "AnnotationSetsReactome", function(.Object, genesets_table,
                                                     genesets_proteins,
                                                     versions) {

  # Define the required column names
  required_genesets_table <- c(getSharedColNamesAnnotations(), "date_download" )

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

  return(.Object)
})

#######################
##### IntensitiesTable
#' @export
setClass("IntensityDataset", slots = c(intensitiesTable = "data.frame", metadataTable = "data.frame"))

setMethod("initialize", "IntensityDataset", function(.Object, intensitiesTable,
                                                     metadataTable) {

  # Define the required column names
  required_intensity_columns <- c("GroupId", "replicate", "NormalisedIntensity", "Imputed",
                                  "replicateNumber", "condition", "ReplicateCountTot", "numberOfReplicateCount",
                                  "precentageOfReplicates")

  # Check if the input data frame has the required column names
  if (!all(required_intensity_columns %in% colnames(intensitiesTable))) {
    stop("intensitiesTable must have the following column names: ", paste(required_intensity_columns, collapse = ", "))
  }

  # Define the required column names
  required_metadata_columns <- c("GroupId", "GroupLabel", "ProteinIds", "ProteinNames", "ProteinGroup", "GeneNames",
                     "Description","QValue", "GroupLabelType")
  if (!all(required_metadata_columns %in% colnames(metadataTable))) {
    stop("metadataTable must have the following column names: ", paste(required_metadata_columns, collapse = ", "))
  }

  groupIdsIntensity <- unique(intensitiesTable$GroupId)
  groupIdsMetadata <- unique(metadataTable$GroupId)
  if(!all(groupIdsIntensity %in% groupIdsMetadata)){
    stop("Some GroupIds in the intensitiesTable are missing from the metadataTable.")
  }

  # Assign the data frame to the data slot
  .Object@intensitiesTable <- intensitiesTable
  .Object@metadataTable <- metadataTable

  return(.Object)
})


#' IntensityDataset
#' @export
createIntensityDataset <- function(intensitiesTable, metadataTable) {
  return(new("IntensityDataset", intensitiesTable = intensitiesTable, metadataTable = metadataTable))
}

