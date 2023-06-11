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
    if (length(intersect(rowData(comparisonExperiment)$ProteinId, unlist(geneSets)))>100){

      # do enrichment
      if (method == "camera"){
        enrichmentStatistics <- .manualRunCameraEnrichment(comparisonExperiment, geneSets)
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
      print(rowData(comparisonExperiment)$ProteinId[1:10])
      print(unlist(geneSets)[1:10])
    }

  }

  comparisonEnrichment
}

.manualRunCameraEnrichment <- function(comparisonExperiment, geneSets){
  #Part of this code was derived from the EnrichmentBrowser package.


  print("Running CAMERA using native code")


  # Prepare CAMERA Inputs
  expression_data <- assay(comparisonExperiment)
  grp <- colData(comparisonExperiment)$GROUP
  group <- factor(grp)
  f <- "~"
  f <- formula(paste0(f, "group"))
  design <- model.matrix(f)

  # filter for genes in the experiment
  igenes <- intersect(rowData(comparisonExperiment[rownames(expression_data)])$ProteinId, unique(unlist(geneSets)))
  igenes_rownames <- rownames(
    rowData(comparisonExperiment)[rowData(comparisonExperiment)$ProteinId %in% igenes,]
  )

  if(!length(igenes)) stop("Expression dataset (se)", " and ",
                           "gene sets (geneSets) have no gene IDs in common")
  expression_data <- expression_data[igenes_rownames,]
  filtered_geneSets <- lapply(geneSets, function(s) s[s %in% igenes])


  # Call CAMERA
  geneSetsIndex <- limma::ids2indices(filtered_geneSets, igenes)
  cameraEnrichmentResults <- limma::camera(expression_data, geneSetsIndex, design, sort=FALSE)

  # format results
  cameraEnrichmentResults$Pathway <- rownames(cameraEnrichmentResults)

  cameraEnrichmentResults <- cameraEnrichmentResults[, c("Pathway","PValue", "FDR","NGenes")]
  colnames(cameraEnrichmentResults) <- c("gene.set", "pval", "adj.pval", "observed")

  cameraEnrichmentResults
}

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
  mask <- rowData(comparisonExperiment)$ProteinId %in% geneSet[[1]]
  mean(rowData(comparisonExperiment)[mask, "FC"])
}

#' add Average Fold Change statistics to enrichment Results
#' @export addAverageFoldChangeColumn
addAverageFoldChangeColumn <- function(enrichmentResults, comparisonExperiment){
  enrichmentResults$afc <- unlist(lapply(enrichmentResults$items,
                                         function(x){calculateAverageFoldChange(comparisonExperiment, x)}))
  enrichmentResults
}
