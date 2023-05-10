Mass Dynamics enrichment service
================

- <a href="#new-enrichment-workflow" id="toc-new-enrichment-workflow">New
  enrichment workflow</a>
  - <a href="#make-annotations" id="toc-make-annotations">Make
    annotations</a>
  - <a href="#read-dataset-parameters" id="toc-read-dataset-parameters">Read
    dataset parameters</a>

``` r
devtools::install_github("MassDynamics/EnrichmentAnalysisStepR")
```

``` r
library(EnrichmentAnalysisStepR)
utils::packageVersion("EnrichmentAnalysisStepR")
```

    ## [1] '0.0.27'

# New enrichment workflow

## Make annotations

- Currently missing Chinese Hamster

This will build GO, Reactome for Human, Mouse, Yeast and MsigDB
annotations sets for Human as internal data to the package stored in
`inst/extdata`

``` bash
make all BUILD_PACKAGE=True MSIG_VERSION="7.5"
```

List data available in the package

``` r
extdata_path <- system.file("extdata", package = "EnrichmentAnalysisStepR")
data_files <- list.files(extdata_path, recursive = TRUE, pattern = "*rda")

data_files
```

    ##  [1] "10090/GO_mouse_2023-05-10.rda"      "10090/Reactome_2023-05-10.rda"     
    ##  [3] "559292/GO_yeast_2023-05-10.rda"     "559292/Reactome_2023-05-10.rda"    
    ##  [5] "9606/GO_human_2023-05-10.rda"       "9606/msigdb_v7.5_2023-05-10.rda"   
    ##  [7] "9606/msigdb_v7.5_C1_2023-05-10.rda" "9606/msigdb_v7.5_C2_2023-05-10.rda"
    ##  [9] "9606/msigdb_v7.5_C3_2023-05-10.rda" "9606/msigdb_v7.5_C4_2023-05-10.rda"
    ## [11] "9606/msigdb_v7.5_C5_2023-05-10.rda" "9606/msigdb_v7.5_C6_2023-05-10.rda"
    ## [13] "9606/msigdb_v7.5_C7_2023-05-10.rda" "9606/msigdb_v7.5_C8_2023-05-10.rda"
    ## [15] "9606/msigdb_v7.5_H_2023-05-10.rda"  "9606/Reactome_2023-05-10.rda"

Example load gene set.

``` r
load(system.file("extdata","9606/GO_human_2023-05-10.rda", package = "EnrichmentAnalysisStepR"))
head(go_sets$genesets_table)
```

    ## # A tibble: 6 × 7
    ##   annotation_id annotation_name    annotation_subcategory annotation_description
    ##   <chr>         <chr>              <chr>                  <chr>                 
    ## 1 GO:0003723    RNA binding        molecular_function     Binding to an RNA mol…
    ## 2 GO:0046872    metal ion binding  molecular_function     Binding to a metal io…
    ## 3 GO:0005829    cytosol            cellular_component     The part of the cytop…
    ## 4 GO:0002250    adaptive immune r… biological_process     An immune response me…
    ## 5 GO:0005886    plasma membrane    cellular_component     The membrane surround…
    ## 6 GO:0019814    immunoglobulin co… cellular_component     A protein complex tha…
    ## # ℹ 3 more variables: size <int>, linkout <chr>, date_download <date>

## Read dataset parameters

``` r
# example read json file
gene_sets <- list("BiologicalProcess", "MolecularFunction")
condition_comparisons <- list(list(right = "COVID", left = "MP"), 
                              list(right = "COVID", left = "Control"))
datasetParams <- list(dataset_id = "ID1", 
                      version = utils::packageVersion("EnrichmentAnalysisStepR"), 
                      fields = list(experiment_design = data.frame(), 
                                    species_gene_set_selection = list(gene_sets=gene_sets, species = "Human"), 
                                                  condition_comparisons = condition_comparisons))
```

``` r
# create class 
datasetParams <- new("EnrichmentDatasetParameters", dataset_id = datasetParams$dataset_id,
                     version = datasetParams$dataset_id, fields = datasetParams$fields)

getConditionComparisons(datasetParams)
```

    ## # A tibble: 2 × 2
    ##   right left   
    ##   <chr> <chr>  
    ## 1 COVID MP     
    ## 2 COVID Control

``` r
getSpecies(datasetParams)
```

    ## [1] "Human"

``` r
getGeneSets(datasetParams)
```

    ## [1] "BiologicalProcess" "MolecularFunction"
