library(testthat)
library(EnrichmentAnalysisStepR)


test_proteinids_parser <- function(){

  real_uniprot_one_protein_sp <- "sp|FLT3|ACN2_ACAGO"
  real_uniprot_one_protein_tr <- "tr|FLT3|ACN2_ACAGO"

  real_uniprot_two_proteins_sp <- "sp|FLT3;P27748|ACN2_ACAGO"
  real_uniprot_two_proteins_tr <- "tr|FLT3;P27748|ACN2_ACAGO"

  not_uniprot_with_semicolon <- "aa|333;bb|444"
  not_uniprot <- "aa|333|444"


  test_that("Real uniprot IDs with one protein are parsed correctly", {

    outcome_sp <- parsePipeId(real_uniprot_one_protein_sp)
    outcome_tr <- parsePipeId(real_uniprot_one_protein_tr)
    outcome_simple <- parsePipeId("FLT3")
    expect_equal(outcome_sp, "FLT3")
    expect_equal(outcome_tr, "FLT3")
    expect_equal(outcome_simple, "FLT3")

  })

  test_that("Real uniprot IDs with two proteins are parsed correctly", {

    outcome_sp <- parsePipeId(real_uniprot_two_proteins_sp)
    outcome_tr <- parsePipeId(real_uniprot_two_proteins_tr)
    expect_equal(outcome_sp, "FLT3")
    expect_equal(outcome_tr, "FLT3")

  })

  test_that("Non uniprot IDs with with semicolon in the string are parsed correctly", {

    outcome <- parsePipeId(not_uniprot_with_semicolon)
    expect_equal(outcome, "aa|333;bb|444")

  })

  test_that("Non uniprot IDs with pipes are parsed correctly", {

    outcome <- parsePipeId(not_uniprot)
    expect_equal(outcome, not_uniprot)

  })

}


test_proteinids_parser()

