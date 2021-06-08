library(EnrichmentAnalysisStepR)
library(testthat)
library(jsonlite)

acceptance_test<- function(path_to_test_data,
                           path_to_expected_data,
                           tolerance = 10**-3){

  results = enrichment_workflow_step(path_to_test_data,
                                     gmt_folder = path_to_test_data,
                                     output_folder = path_to_test_data,
                                     by = "Protein",
                                     method = "camera")

  current = read_json(file.path(path_to_test_data, "C.P.Reactome.json"))
  expected = read_json(file.path(path_to_test_data, "expected.C.P.Reactome.json"))

  test_that("approximately equal", {
    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })

}


path = "../../test_data/"
acceptance_test(path,path)
