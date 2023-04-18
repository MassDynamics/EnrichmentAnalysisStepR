library(EnrichmentAnalysisStepR)

test_datasets <- function() {
  datasetParametersPath <- "../test_data/datasets/example_dataset_parameters.json"
  datasetParameters <- jsonlite::fromJSON(datasetParametersPath)

  test_getConditionComparisons <- (function() {
    test_that("returns the comparisons required for the dataset", {
      conditionComparisons <- EnrichmentAnalysisStepR::datasets.getConditionComparisons(
        datasetParameters
      )

      expectedConditionComparisons <- data.frame(
        left = c("S1", "S2", "S1"), 
        right = c("zero", "zero", "S2")
      )

      expect_equal(conditionComparisons, expectedConditionComparisons)
    })
  })()
  
  test_getGeneSetsInfo <- (function() {
    test_that("returns the comparisons required for the dataset", {
      geneSetsInfo <- EnrichmentAnalysisStepR::datasets.getGeneSetsInfo(
        datasetParameters
      )

      expectedGeneSetsInfo <- list(
        species = "Human", 
        gene_set_names = c("BiologicalProcess", "MolecularFunction")
      )

      expect_equal(geneSetsInfo, expectedGeneSetsInfo)
    })
  })()
}

test_datasets()