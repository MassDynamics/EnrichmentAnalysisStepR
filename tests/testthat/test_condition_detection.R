library(EnrichmentAnalysisStepR)
library(testthat)
library(jsonlite)

test_condition_detection <- function(){

  
  test_that("simple condition match", {
    
    conditions = c("Parental","SKBR3")
    comparison = "Parental SKBR3"
    
    second_condition = get_condition_string(conditions,comparison, 2)
    first_condition = get_condition_string(conditions,comparison, 1)
    
    expect_true(second_condition == "SKBR3") 
    expect_true(first_condition == "Parental")
  })
  
  test_that("condition match with special characters",{
    
    conditions = c("AZD8931-resistant SKBR3-AZDRc","Parental SKBR3")
    comparison = "AZD8931-resistant SKBR3-AZDRc - Parental SKBR3"
    
    second_condition = get_condition_string(conditions,comparison, 2)
    first_condition = get_condition_string(conditions,comparison, 1)
    
    expect_true(second_condition == "Parental SKBR3") 
    expect_true(first_condition == "AZD8931-resistant SKBR3-AZDRc")
  })
}

test_condition_detection()