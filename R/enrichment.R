#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#' @param comparison A summarized experiment object containing only two experimental groups
#' @param gmtFolder a folder containing gene set libraries which need to be referenced
#' @param method the method enrichment browser should use for set based enrichment
#' @return enrichment An enrichment table provided by EnrichmentBrowserPackage
#' @export enrichComparisonExperiment
enrichComparisonExperiment <- function(comparisonExperiment, gmtFolder, method, perm = perm){

  comparisonEnrichment <- list(libraryStatistics = list(),
                               up.condition = metadata(comparisonExperiment)$up.condition,
                               down.condition = metadata(comparisonExperiment)$down.condition)

  for (library in list.files(gmtFolder, pattern = "\\.gmt$")){
    gmt.file <- file.path(gmtFolder, library)

    geneSets <- getGenesets(gmt.file)

    print(gmt.file)
    print(paste("Sets in Library: ", length(geneSets)))
    print(paste("Genes in Library: ", length(unique(unlist(geneSets)))))
    print(paste("Genes in Library and dataset: ", length(intersect(rownames(comparisonExperiment),unlist(geneSets)))))

    # check for common ids
    if (length(intersect(rownames(comparisonExperiment),unlist(geneSets)))>100){

      # do enrichment
      if (method == "camera"){
        enrichmentStatistics <- manualRunCameraEnrichment(comparisonExperiment, geneSets)
      } else {
        sbeaResults <- sbea(method = method, se = comparisonExperiment,
                            geneSets = gmt.file, perm = perm, alpha = 1)
        enrichmentStatistics = as.data.frame(gsRanking(sbeaResults))
        colnames(enrichmentStatistics) <- tolower(colnames(enrichmentStatistics))
      }

      print(colnames(enrichmentStatistics))
      # pval and gene.set must be present!
      for (col in c("gene.set", "pval")){
        stopifnot(col %in% colnames(enrichmentStatistics))
      }

      # write enrichmentStatistics
      lib = str_split(library,".gmt")[[1]][1]
      label = str_c(lib,"_enrichment.csv")

      enrichmentStatistics <- addGeneSetsColumn(enrichmentStatistics, geneSets)
      enrichmentStatistics <- addAverageFoldChangeColumn(enrichmentStatistics, comparisonExperiment)

      # add adjusted p value
      enrichmentStatistics$adj.pval = p.adjust(enrichmentStatistics$pval, method = "BH")

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


#' Wrapper for CAMERA
#' @param comparisonExperiment A summarized experiment object containing only two experimental groups
#' @param geneSets gene sets list to run CAMERA
#' @export manualRunCameraEnrichment
manualRunCameraEnrichment <- function(comparisonExperiment, geneSets){
  #Part of this code was derived from the EnrichmentBrowser package.

  print("Running CAMERA using native code")

  # Prepare CAMERA Inputs
  expression_data <- assay(comparisonExperiment)
  grp <- colData(comparisonExperiment)$GROUP
  group <- factor(grp)
  f <- formula("~group")
  design <- model.matrix(f)

  # filter for genes in the experiment
  igenes <- intersect(rownames(expression_data), unique(unlist(geneSets)))
  if(!length(igenes)) stop("Expression dataset (se)", " and ",
                           "gene sets (geneSets) have no gene IDs in common")
  expression_data <- expression_data[igenes,]
  filtered_geneSets <- lapply(geneSets, function(s) s[s %in% igenes])


  # Call CAMERA
  geneSetsIndex <- limma::ids2indices(filtered_geneSets, rownames(expression_data))
  cameraEnrichmentResults <- limma::camera(expression_data, geneSetsIndex, design, sort=FALSE)

  # format results
  cameraEnrichmentResults$Pathway <- rownames(cameraEnrichmentResults)

  cameraEnrichmentResults <- cameraEnrichmentResults[, c("Pathway","PValue", "FDR","NGenes")]
  colnames(cameraEnrichmentResults) <- c("gene.set", "pval", "adj.pval", "observed")

  cameraEnrichmentResults
}


#' Wrapper for CAMERA
#' @param comparisonExperiment A summarized experiment object containing only two experimental groups
#' @param geneSets gene sets list to run CAMERA
#' @preRanked logical. TRUE to use the pre-ranked version of `camera`.
#' @export runCamera
runCamera <- function(comparisonExperiment, geneSets, preRanked= FALSE){
  #Part of this code was derived from the EnrichmentBrowser package.

  print("Running CAMERA using native code")

  # Prepare CAMERA Inputs
  expression_data <- assay(comparisonExperiment)
  grp <- colData(comparisonExperiment)$GROUP
  group <- factor(grp)
  f <- formula("~group")
  design <- model.matrix(f)

  # filter for genes in the experiment
  igenes <- intersect(rownames(expression_data), unique(unlist(geneSets)))
  if(!length(igenes)) stop("Expression dataset (se)", " and ",
                           "gene sets (geneSets) have no gene IDs in common")

  geneSetsIndex <- limma::ids2indices(geneSets, rownames(expression_data))

  # Call CAMERA
  if(preRanked){
    fit <- eBayes(lmFit(expression_data, design))
    cameraEnrichmentResults <- cameraPR(fit$t[,2], geneSetsIndex, sort = FALSE)
  } else {
    cameraEnrichmentResults <- limma::camera(expression_data, geneSetsIndex, design, sort=FALSE, trend.var = TRUE)
  }

  # format results
  cameraEnrichmentResults$Pathway <- rownames(cameraEnrichmentResults)

  cameraEnrichmentResults <- cameraEnrichmentResults[, c("Pathway","PValue", "FDR","NGenes","Direction")]
  colnames(cameraEnrichmentResults) <- c("gene.set", "pval", "adj.pval", "observed","direction")

  cameraEnrichmentResults
}


# utility for adding gene set column name
#' @export addGeneSetsColumn
addGeneSetsColumn <- function(result_table, geneSets){
  result_table$items <- geneSets[match(result_table$gene.set, names(geneSets))]
  result_table
}

# utility for converting list of char vectors to table
#' @export geneSetstoDataFrame
geneSetstoDataFrame <- function(geneSets){
  df <- as.data.frame(cbind(names(geneSets),geneSets))
  colnames(df) <- c("gene.set","items")
  df$gene.set <- sapply(df$gene.set, toString)
  df
}

# get's average fold change using a gene set from the limma fold change statistics
# stored in the rowData for the comparison experiments
#' @export calculateAverageFoldChange
calculateAverageFoldChange <- function(comparisonExperiment, geneSet){
  mean(rowData(comparisonExperiment)[rownames(comparisonExperiment) %in% geneSet,"FC"])
}

#' add Average Fold Change statistics to enrichment Results
#' @export addAverageFoldChangeColumn
addAverageFoldChangeColumn <- function(enrichmentResults, comparisonExperiment){
  enrichmentResults$afc <- unlist(lapply(enrichmentResults$items,
                                         function(x){calculateAverageFoldChange(comparisonExperiment, x)}))
  enrichmentResults
}
