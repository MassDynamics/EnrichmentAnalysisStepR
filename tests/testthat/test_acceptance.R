library(EnrichmentAnalysisStepR)
library(testthat)
library(jsonlite)

acceptance_test<- function(path_to_test_data,
                           path_to_expected_data,
                           tolerance = 10**-3){

  results = runEnrichmentWorkflowStep(path_to_test_data,
                                     gmtFolder = path_to_test_data,
                                     outputFolder = path_to_test_data,
                                     by = "Protein",
                                     method = "camera")

  current = read_json(file.path(path_to_test_data, "C___P_Reactome.json"))
  expected = read_json(file.path(path_to_test_data, "expected_CPReactome.json"))

  test_that("approximately equal", {
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })

}


path = "../test_data/"
acceptance_test(path, path)
