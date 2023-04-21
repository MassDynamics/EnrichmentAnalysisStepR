#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#' @param comparison A summarized experiment object containing only two experimental groups
#' @param gmtFolder a folder containing gene set libraries which need to be referenced
#' @param method the method enrichment browser should use for set based enrichment
#' @return enrichment An enrichment table provided by EnrichmentBrowserPackage
#' @export runComparisonExperiment
runComparisonExperiment <- function(intensitiesDf, metadataDf,
                                    conditionComparison,
                                    geneSetsInfo,
                                    outputFormat = "tsv",
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
                                                      metadataDf,
                                                      conditionComparison,
                                                      conditionCol,
                                                      replicateCol,
                                                      featureIdCol,
                                                      intensityCol,
                                                      imputedCol)

  # add stop if species is not in the available data sets
  species <- geneSetsInfo$species
  libraries <- geneSetsInfo$gene_set_names

  comparisonEnrichment <- list(libraryStatistics = list(),
                               up.condition = conditionComparison$rightCondition,
                               down.condition = conditionComparison$leftCondition)

  for (library in libraries){
    # make it a package
    geneSets <- GENE_SETS_LIST[[species]][[library]]

    print(glue("Species: {species}; Gene set: {library}"))
    print(glue("Sets in Library: {length(geneSets)}"))
    print(glue("Genes in Library: {length(unique(unlist(geneSets)))}"))
    print(glue("Genes in Library and dataset: {length(intersect(rownames(comparisonExperiment),unlist(geneSets)))}"))

    # check for common ids
    if (length(intersect(rownames(comparisonExperiment),unlist(geneSets)))>100){

      # do enrichment
      enrichmentStatistics <- runCamera(comparisonExperiment, geneSets, method = method)

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

      comparisonEnrichment$libraryStatistics[[lib]] = enrichmentStatistics

    } else {
      print("Skipping enrichment, not enough genes in common.")
      print("Heres some example ids")
      print(rownames(comparisonExperiment)[1:10])
      print(unlist(geneSets)[1:10])
    }

  }

  comparisonEnrichment
}



#' Orchestrates the enrichment workflow step for MD input/output data
#' @export runEnrichmentWorkflowStep
runEnrichmentWorkflowStep <- function(uploadFolder,
                                      gmtFolder,
                                      outputFolder,
                                      by = "Protein",
                                      method = "camera",
                                      perm){

  if (file.exists(file.path(uploadFolder,"DiscoveryQuant.RData"))){
    print("Found Discovery Data Export. Loading data")
    load(file.path(uploadFolder,"DiscoveryQuant.RData"))
  } else {
    print("No Discovery Data Export Found. Using Object Export")
    mdQuantExports <- loadMDQuantOutput(uploadFolder)
    comparisonExperiments <- listcomparisonExperimentsList(mdQuantExports)
  }

  enrichmentResults <- enrichComparisons(
    comparisonExperiments,
    gmtFolder
  )

  enrichmentResults <- annotateEnrichmentResults(enrichmentResults,
                                                 gmtFolder)
  writeEnrichmentJSON(enrichmentResults, outputFolder)

  # done!
  enrichmentResults
}

#' Orchestrates the enrichment workflow step for a comparisonExperiment List
#' @export enrichComparisons
enrichComparisons <- function(comparisonExperimentsList,
                              gmtFolder,
                              method = "camera",
                              perms = 100){

  results = list()
  for (comparison in names(comparisonExperimentsList)){
    print(paste0("Running enrichment for comparison: ", comparison))
    comparisonExperiment = comparisonExperimentsList[[comparison]]
    results[[comparison]] = enrichComparisonExperiment(comparisonExperiment, gmtFolder, method, perms)
  }

  results
}
