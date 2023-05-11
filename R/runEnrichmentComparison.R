#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#' @param comparison A summarized experiment object containing only two experimental groups
#' @param gmtFolder a folder containing gene set libraries which need to be referenced
#' @param method the method enrichment browser should use for set based enrichment
#' @return enrichment An enrichment table provided by EnrichmentBrowserPackage
#' @export
runEnrichmentOneAnnotation <- function(intensityDataset,
                                    conditionComparisons,
                                    geneSetsList,
                                    database,
                                    method = "camera",
                                    conditionCol = "condition",
                                    replicateCol = "replicate",
                                    featureIdCol = "GroupId",
                                    intensityCol = "NormalisedIntensity",
                                    imputedCol = "Imputed"){


  nComparisons <- nrow(conditionComparisons)
  info(EnrichmentLogger(), glue("Number of pairwise comparisons: {nComparisons}"))

  enrichmentComparisonList <- list()
  for(comparisonIdx in 1:nrow(conditionComparisons)){
    conditionComparison <- conditionComparisons[comparisonIdx,]
    info(EnrichmentLogger(), glue("Running {comparisonIdx} comparison: left {conditionComparison$left}; right {conditionComparison$right}"))

    comparisonExperiment <- prepareComparisonExperiment(intensityDataset,
                                                        conditionComparison,
                                                        conditionCol,
                                                        replicateCol,
                                                        featureIdCol,
                                                        intensityCol,
                                                        imputedCol)

    info(EnrichmentLogger(), "comparisonExperiment ready for enrichment analysis.")
    info(EnrichmentLogger(), glue("Genes in Library and dataset: {length(intersect(rownames(comparisonExperiment),unlist(geneSetsList)))}"))

    # check for common ids
    if (length(intersect(rownames(comparisonExperiment),unlist(geneSetsList)))>100){

      # do enrichment
      enrichmentStatistics <- runGeneSetTest(comparisonExperiment, geneSetsList, method)

      # pval and gene.set must be present!
      for (col in c("gene.set", "pval")){
        stopifnot(col %in% colnames(enrichmentStatistics))
      }

      # add these - annotate using the other file and
      enrichmentStatistics <- addGeneSetsColumn(enrichmentStatistics, geneSetsList)
      enrichmentStatistics <- addAverageFoldChangeColumn(enrichmentStatistics, comparisonExperiment)
      enrichmentStatistics$database <- database
      enrichmentStatistics$leftCondition <- conditionComparison$left
      enrichmentStatistics$rightCondition <- conditionComparison$right

      stopifnot("afc" %in% colnames(enrichmentStatistics))
      stopifnot("items" %in% colnames(enrichmentStatistics))

      enrichmentComparisonList[[comparisonIdx]] <- enrichmentStatistics

    } else {
      print("Skipping enrichment, not enough genes in common.")
      print("Heres some example ids")
      print(rownames(comparisonExperiment)[1:10])
      print(unlist(geneSetsList)[1:10])
    }
  }
  comparisonEnrichment <- bind_rows(enrichmentComparisonList)
  return(comparisonEnrichment)
}


#' Automates Enrichment Analyses for one dataaset
#' @param intensitiesDf List of three list each containing a pair of names and tables for: protein metadata,
#' intensities, condition comparison
#' @param metadataDf Features metadata
#' @param datasetParams list of metadata passed by the user
#' @import glue

#' @export
runEnrichmentForAnnotationSet <- function(intensitiesDf, metadataDf, datasetParams){

  intensityDataset <- createIntensityDataset(intensitiesTable = intensitiesDf, metadataTable = metadataDf)
  datasetParams <- createEnrichmentDatasetParameters(dataset_id = datasetParams$dataset_id,
                                                     version = datasetParams$version,
                                                     fields = datasetParams$fields)

  species <- getSpecies(datasetParams)
  conditionComparisons <- getConditionComparisons(datasetParams)
  geneSets <- getGeneSets(datasetParams)
  annotationDate <- getAnnotationDate()
  info(EnrichmentLogger(), glue("Download date of annotations currently used: {annotationDate}"))

  if(is.null(species)){
    stop("Species is missing") # TODO : when other
  }

  # for each gene set - run through all the pairiwse comparisons
  resultsList <- list()
  for (geneSet in geneSets){
    currentSet <- getGeneSetList(species = species, geneSet = geneSet, version = annotationDate)
    currentSetList <- df_items_to_list(currentSet@genesets_proteins)

    info(EnrichmentLogger(), glue("Species: {species}; Gene set: {geneSet}"))
    info(EnrichmentLogger(), glue("Sets in Library: {length(currentSetList)}"))
    info(EnrichmentLogger(), glue("Genes in Library: {length(unique(unlist(currentSetList)))}"))

    results <- runEnrichmentOneAnnotation(intensityDataset = intensityDataset,
                                          geneSetsList = currentSetList,
                                          database = geneSet,
                                          conditionComparison = conditionComparisons)

    results <- results |> left_join(currentSet@genesets_table, by = c("gene.set" = "annotation_id"))
    resultsList[[geneSet]] <- results
  }

  results <- bind_rows(resultsList)
  results <- results |> dplyr::rename(annotation_id = gene.set)
  return(results)

}

#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#' @param species species, e.g. Human, Mouse
#' @param geneSet name of gene set: Reactome, MolecularProcess etc..
#' @import glue
#' @export
getGeneSetList <- function(species, geneSet, version){
  taxonomy <- species_to_taxonomy(species)

  if(species == "Human"){
    setFile <- map_sets_to_file_human(geneSet)
  }

  if(species == "Mouse"){
    setFile <- map_sets_to_file_mouse(geneSet)
  }

  if(species == "Yeast"){
    setFile <- map_sets_to_file_yeast(geneSet)
  }
  setFilePath <- glue("{version}/{taxonomy}/{setFile}")

  assign("currentSetName", load(system.file("extdata", setFilePath, package = "EnrichmentAnalysisStepR")))
  return(get(currentSetName))
}


#' Get date of annotation currently used
#' @export
getAnnotationDate <- function(){
  version <- Sys.getenv("ANNNOTATION_DATE")
  # If the variable is not set, an empty string will be returned
  if (version == "") {
    stop("Environment variable ANNNOTATION_DATE is not set.")
  }
  return(version)
}

#datasetParams <- jsonlite::read_json("development/trigger_enrich_manually/enrichment_dataset_example.json")
#enrichParams <- new("EnrichmentDatasetParameters", dataset_id=dataset_params$dataset_id,
#                    version = dataset_params$version, fields = dataset_params$fields)
