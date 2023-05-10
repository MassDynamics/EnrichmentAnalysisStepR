

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
