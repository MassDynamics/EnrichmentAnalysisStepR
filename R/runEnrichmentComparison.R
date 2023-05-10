#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#' @param comparison A summarized experiment object containing only two experimental groups
#' @param gmtFolder a folder containing gene set libraries which need to be referenced
#' @param method the method enrichment browser should use for set based enrichment
#' @return enrichment An enrichment table provided by EnrichmentBrowserPackage
#' @export runEnrichmentComparison
runEnrichmentComparison <- function(intensitiesDf,
                                    conditionComparison, # TODO: accept multiple comparisons - multiple rows
                                    geneSetsInfo,
                                    method = "camera",
                                    conditionCol = "condition",
                                    replicateCol = "replicate",
                                    featureIdCol = "GroupId",
                                    intensityCol = "NormalisedIntensity",
                                    imputedCol = "Imputed"){

  if(nrow(conditionComparison) > 1){
    warning("conditionComparison has more than one row and only the first pairwise comparison will be considered.")
  }

  comparisonExperiment <- prepareComparisonExperiment(intensitiesDf,
                                                      conditionComparison,
                                                      conditionCol,
                                                      replicateCol,
                                                      featureIdCol,
                                                      intensityCol,
                                                      imputedCol)

  comparisonEnrichment <- list(libraryStatistics = list(),
                               up.condition = conditionComparison$rightCondition,
                               down.condition = conditionComparison$leftCondition)

  for (library in libraries){
    # MAKE IT A PACKAGE
    warning("I'M USING A FAKE LIBRARY - FIX THIS IN RUNNER.R")
    geneSets <- GENE_SETS_LIST[[species]][[library]]

    print(glue("Species: {species}; Gene set: {library}"))
    print(glue("Sets in Library: {length(geneSets)}"))
    print(glue("Genes in Library: {length(unique(unlist(geneSets)))}"))
    print(glue("Genes in Library and dataset: {length(intersect(rownames(comparisonExperiment),unlist(geneSets)))}"))

    # check for common ids
    if (length(intersect(rownames(comparisonExperiment),unlist(geneSets)))>100){

      # do enrichment
      enrichmentStatistics <- runGeneSetTest(comparisonExperiment, geneSets, method)

      print(colnames(enrichmentStatistics))
      # pval and gene.set must be present!
      for (col in c("gene.set", "pval")){
        stopifnot(col %in% colnames(enrichmentStatistics))
      }

      # add these - annotate using the other file and
      enrichmentStatistics <- addGeneSetsColumn(enrichmentStatistics, geneSets)
      enrichmentStatistics <- addAverageFoldChangeColumn(enrichmentStatistics, comparisonExperiment)

      print(colnames(enrichmentStatistics))
      stopifnot("afc" %in% colnames(enrichmentStatistics))
      stopifnot("items" %in% colnames(enrichmentStatistics))

      comparisonEnrichment$libraryStatistics[[library]] = enrichmentStatistics

    } else {
      print("Skipping enrichment, not enough genes in common.")
      print("Heres some example ids")
      print(rownames(comparisonExperiment)[1:10])
      print(unlist(geneSets)[1:10])
    }

  }

  comparisonEnrichment <- annotateComparison(comparisonEnrichment)

  comparisonEnrichment
}


#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#' @param datasetsTables List of three list each containing a pair of names and tables for: protein metadata,
#' intensities, condition comparison
#' @param metadata list of metadata passed by the user

#' @export runEnrichmentForAnnotationSet
runEnrichmentForAnnotationSet <- function(datasetParams){

  datasetParams <- new()

  if(!is.list(geneSetsInfo) & !(names(geneSetsInfo) %in% c("species", "")))
  # add stop if species is not in the available data sets
  species <- geneSetsInfo$species
  libraries <- geneSetsInfo$gene_set_names
  }



dataset_params <- jsonlite::read_json("development/trigger_enrich_manually/enrichment_dataset_example.json")
#enrichParams <- new("EnrichmentDatasetParameters", dataset_id=dataset_params$dataset_id,
#                    version = dataset_params$version, fields = dataset_params$fields)
