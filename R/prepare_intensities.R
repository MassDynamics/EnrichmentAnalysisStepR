#' Prepare comparison experiment
#'
#' This function prepares a comparison experiment by organizing input data into a standardized format.
#'
#' @param intensitiesDf A data.frame containing intensity values for features.
#' @param metadataDf A data.frame containing metadata for each feature.
#' @param conditionComparison list containing levels of conditions to be compared.
#' @param conditionCol A character string specifying the column name in intensitiesDf that contains the condition information.
#' @param replicateCol A character string specifying the column name in intensitiesDf that contains the replicate information.
#' @param featureIdCol A character string specifying the column name in intensitiesDf and metadataDf that contains the feature identifier information.
#' @param intensityCol A character string specifying the column name in intensitiesDf that contains the intensities.
#' @param imputedCol A character string specifying the column name in intensitiesDf that contains the imputed value information.
#' @param replicateConditionMapping A list specifying the mapping between replicates and conditions.
#'
#' @return a SummarizedExperiment object containing the intensities, column and row data to perform gene set testing.

#' @importFrom dplyr left_join
#' @export
#'
# the input could be simplified if we create a fixed class for intensitiesDf
prepareComparisonExperiment <- function(intensitiesDf,
                                        metadataDf,
                                        conditionComparison,
                                        conditionCol,
                                        replicateCol,
                                        featureIdCol,
                                        intensityCol,
                                        imputedCol){

  leftCondition <- conditionComparison$left[1]
  rightCondition <- conditionComparison$right[1]

  # Replicates and conditions
  replicateConditionMapping <- unique(intensitiesDf[,c(conditionCol,replicateCol)])
  uniqueConditions <- unique(replicateConditionMapping[,conditionCol])
  if(!(leftCondition %in% uniqueConditions) | !(rightCondition %in% uniqueConditions)){
    stop("The up and down conditions provided are not present in the samples data.")
  }

  filtered_data <- filterByMissingProportion(longIntensityDT = intensitiesDf,
                                             conditionColname = conditionCol,
                                             runIdColname = replicateCol,
                                             idCol = featureIdCol,
                                             imputedColname = imputedCol,
                                             minPropAvail = 0.5)

  filtered_data <- filtered_data %>% left_join(metadataDf)
  proteinIntensitiesWide <- makeWideIntensitiesDf(filtered_data, groupByCol = "ProteinIds",
                                                  intensityCol = intensityCol, replicateCol = replicateCol)

  comparisonExperiment <- createComparisonSummarisedExperiment(proteinIntensitiesWide,
                                                               metadataDf=metadataDf,
                                                               idCol = "ProteinIds",
                                                               conditionColname = conditionCol,
                                                               leftCondition = leftCondition,
                                                               rightCondition = rightCondition,
                                                               replicateConditionMapping=replicateConditionMapping,
                                                               replicateColname = replicateCol)
  return(comparisonExperiment)
}



#' Convert protein intensity data from long to wide format
#'
#' Given a data frame of protein intensity values in long format, this function
#' converts it to wide format, with one row per protein and one column per
#' replicate. The wide format is useful for downstream analysis and visualization.
#'
#' @param proteinIntensities A data frame containing protein intensity in long format as imported from "Protein_intensity.tsv"
#' @param groupByCol The name of the column in \code{proteinIntensities} that specifies the grouping variable
#'                   to be used for creating the wide format. Default is "ProteinIds".
#' @param intensityCol The name of the column in \code{proteinIntensities} that contains the protein intensity values.
#'                     Default is "NormalisedIntensity".
#' @param replicateCol The name of the column in \code{proteinIntensities} that specifies the replicates.
#'                     Default is "replicate".
#'
#' @return A data frame containing the protein intensity values in wide format
#'
#' @import data.table
#' @importFrom dplyr group_by
#' @importFrom tidyr spread
#'
#' @export makeWideIntensitiesDf
makeWideIntensitiesDf <- function(proteinIntensities, groupByCol = "ProteinIds",
                                  intensityCol = "NormalisedIntensity", replicateCol = "replicate"){

  proteinIntensities <- data.frame(proteinIntensities)
  proteinIntensities <- proteinIntensities[,c(groupByCol, intensityCol, replicateCol)]
  proteinIntensities <- proteinIntensities %>% dplyr::group_by(!!as.symbol(groupByCol), !!as.symbol(replicateCol))
  proteinIntensitiesWide <- proteinIntensities %>%
    dplyr::group_by(!!as.symbol(groupByCol)) %>%
    tidyr::spread(!!as.symbol(replicateCol),  !!as.symbol(intensityCol))

  colnames(proteinIntensitiesWide)[1] = groupByCol
  return(proteinIntensitiesWide)
}


#' Create a SummarizedExperiment object from wide-format protein intensity data
#'
#' Given a data frame of protein intensity values in wide format, this function
#' converts it to a SummarizedExperiment object, which is a Bioconductor class
#' for storing and manipulating genomic data. The resulting object can be used
#' for downstream analysis and visualization.
#'
#' @param wideIntensities A data frame containing protein intensity values in wide format as created from `makeWideIntensitiesDf`
#'
#' @return A SummarizedExperiment object containing the protein intensity values
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export createComparisonSummarisedExperiment
createComparisonSummarisedExperiment <- function(wideIntensities,
                                                 metadataDf,
                                                 leftCondition,
                                                 rightCondition,
                                                 replicateConditionMapping,
                                                 idCol,
                                                 conditionColname,
                                                 replicateColname){

  replicateConditionMapping <- data.frame(replicateConditionMapping)
  replicateMapping <- replicateConditionMapping %>% dplyr::rename(replicate = !!as.symbol(replicateColname),
                                                         condition = !!(as.symbol(conditionColname)))

  replicateMapping <- filterReplicateByConditions(replicateConditionMapping = replicateConditionMapping,
                                                  leftCondition = leftCondition,
                                                  rightCondition = rightCondition)

  replicateMapping <- createColumnData(condition1 = leftCondition, condition2 = rightCondition,
                                       conditionRunIdMapping = replicateMapping, conditionCol = conditionColname)

  wideIntensities <- filterIntensitiesByConditions(wideIntensities = wideIntensities,
                                                   idCol = idCol,
                                                   replicateConditionMapping = replicateMapping)


  intAssay <- createIntensitiesMatrix(wideIntensities, idCol, replicateMapping)

  featureInformation <- prepareRowData(metadataDf, intAssay, idCol, replicateMapping, leftCondition, rightCondition)

  intensitiesSE <- SummarizedExperiment(assays = intAssay, colData = replicateMapping, rowData = featureInformation)

  return(intensitiesSE)
}


#' @export
filterReplicateByConditions <- function(replicateConditionMapping, leftCondition, rightCondition){
  replicateConditionMapping <- replicateConditionMapping[replicateConditionMapping$condition %in% c(leftCondition, rightCondition),]
  if(nrow(replicateConditionMapping) <= 5){
    stop(glue("Not enough samples available. N: {nrow(replicateConditionMapping)}"))
  }

  nLevels <- length(unique(replicateConditionMapping$condition))
  if(nLevels != 2){
    stop(glue("There should only be two levels in the condition under study. N levels detected: {nLevels}"))
  }

  return(replicateConditionMapping)
}


#' @export
filterIntensitiesByConditions <- function(wideIntensities, idCol, replicateConditionMapping){
  requiredColumns = c(idCol,
                      replicateConditionMapping$replicate)
  stopifnot(length(requiredColumns) > 2)

  if(all(!(requiredColumns %in% colnames(wideIntensities)))){
    missing_columns <- paste(requiredColumns[!(requiredColumns %in% colnames(wideIntensities))],collapse=",")
    stop(glue("{missing_columns} columns are missing from intensities table."))
  }

  wideIntensities <- wideIntensities[, requiredColumns]
  return(wideIntensities)
}


#' @export
createIntensitiesMatrix <- function(wideIntensities, idCol, replicateConditionMapping){
  intAssay <- wideIntensities
  intAssay <- intAssay[, -which(colnames(intAssay) == idCol)]
  intAssay <- as.matrix(intAssay)
  rownames(intAssay) <- wideIntensities[,idCol, drop=TRUE]

  orderColnames <- match(colnames(intAssay), replicateConditionMapping$replicate)
  intAssay <- intAssay[,orderColnames]
  stopifnot(colnames(intAssay) == replicateConditionMapping$replicate)
  return(intAssay)
}


#' @importFrom dplyr left_join
#' @export
prepareRowData <- function(metadataDf, intAssay, idCol, replicateMapping, leftCondition, rightCondition){
  rowDataStats <- metadataDf[metadataDf[,idCol] %in% rownames(intAssay),]
  matchOrderIds <- match(rownames(intAssay), rowDataStats[,idCol, drop=TRUE])
  rowDataStats <- rowDataStats[matchOrderIds,]
  stopifnot(rowDataStats[,idCol, drop=TRUE] == rownames(intAssay))

  fcsDf <- computeComparisonFoldChanges(replicateMapping, intAssay, leftCondition, rightCondition, idCol)
  rowDataStatsFCs <- rowDataStats %>% left_join(fcsDf)
  stopifnot(nrow(rowDataStatsFCs) == nrow(rowDataStats))

  return(rowDataStatsFCs)
}


#' @export
computeComparisonFoldChanges <- function(replicateMapping, intAssay, leftCondition, rightCondition, idCol){
  leftCondReplicates <- replicateMapping$replicate[replicateMapping$condition == leftCondition]
  rightCondReplicates <- replicateMapping$replicate[replicateMapping$condition == rightCondition]
  leftCondMeans <- rowMeans(intAssay[,leftCondReplicates])
  rightCondMeans <- rowMeans(intAssay[,rightCondReplicates])
  fcs <- rightCondMeans - leftCondMeans
  fcsDf <- data.frame(FC = fcs, ids = names(fcs))
  colnames(fcsDf)[2] <- idCol
  return(fcsDf)
}
