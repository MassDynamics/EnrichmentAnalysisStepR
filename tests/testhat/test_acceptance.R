library(EnrichmentAnalysisStepR)
library(testthat)

acceptance_test<- function(current, expected,
                           tolerance = 10**-3){

  test_that("columns are in common", {

    result = min(colnames(expected) == colnames(current))
    expect_true(as.logical(result))
  })

  test_that("rows are in common", {

    expected_majority = unique(expected$`majority protein ids`)
    current_majority = unique(current$`majority protein ids`)
    result = length(intersect(expected_majority,current_majority))
    expect_true(result == 3261)
  })

  test_that("approximately equal", {

    approx_same = all.equal(expected, current, tolerance = tolerance)
    expect_true(approx_same) #tolerate small differences
  })

}

output_folder = "../../acceptance_test_data/LFQ/iPRG2015/transform"

# Run Code
tmp =  lfq_transformer(mq_folder = "../../acceptance_test_data/LFQ/iPRG2015/",
                       output_folder = output_folder,
                       imputeStDev=0.3,
                       imputePosition=1.8)
rm(tmp)

path_to_current = file.path(output_folder,"proteinGroups_quant.txt")
current <- fread(path_to_current, sep = "\t", stringsAsFactors = FALSE)
setkey(current, NULL)
path_to_expected = "../../acceptance_test_data/LFQ/iPRG2015/expected_outputs/proteinGroups_quant.txt"
expected <- fread(path_to_expected, sep = "\t", stringsAsFactors = FALSE)
setkey(expected, NULL)

acceptance_test(current, expected)
