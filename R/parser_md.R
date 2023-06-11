
#' Loads the protein quantitation output from an md workflow step
#'
#' @param uploadFolder a folder containing mass dynamics quantitation output
#' @return a named list of mass dyanmics quantitation objects.
#' @export loadMDQuantOutput
loadMDQuantOutput <- function(uploadFolder){

  # Get Experiment Data
  proteinVizPath = file.path(uploadFolder,"protein_viz.json")
  proteinViz = read_json(proteinVizPath, simplifyVector = TRUE)

  proteinIntensitiesLongPath = file.path(uploadFolder,"proteinGroups_quant_intensities.txt")
  proteinIntensitiesLongV2Path = file.path(uploadFolder,"Protein_intensity.tsv")

  if (file.exists(proteinIntensitiesLongPath)) {
    proteinIntensitiesLong = read.csv(proteinIntensitiesLongPath, sep ="\t")
  } else if (file.exists(proteinIntensitiesLongV2Path)) {
    proteinIntensitiesLong = read.csv(proteinIntensitiesLongV2Path, sep ="\t")
    proteinIntensitiesLong = convertV2ProteinIntensties(proteinIntensitiesLong)
  } else {
    proteinIntensitiesLongPath = file.path(uploadFolder,"Protein_Intensity.txt")
    proteinIntensitiesLong = read.csv(proteinIntensitiesLongPath, sep ="\t")
    proteinIntensitiesLong = convertDiscoveryProteinIntensties(proteinIntensitiesLong)
  }

  conditions = as.character(unique(proteinIntensitiesLong$condition))



  mdQuantExports = list(proteinViz = proteinViz,
                        proteinIntensitiesLong = proteinIntensitiesLong,
                        conditions = conditions,
                        comparisons = proteinViz$conditionComparison)

  mdQuantExports
}


#' list all the binary experiments comparisons used for enrichment
#'
#' @param mdQuantExports a named list of mass dyanmics quantitation objects produced by loadMDQuantOutput
#' @return a list of summarized experiments containing protein de expression intensities and statistics
#' @export listcomparisonExperimentsList
listcomparisonExperimentsList <- function(mdQuantExports, conditionComparisonMapping){


  conditions <- mdQuantExports$conditions

  comparisonExperimentsList = list()

  # work out which intensities match which conditions
  conditionComparisonMapping <- getConditionComparisonMapping(mdQuantExports$proteinViz)
  conditionRunIdMapping <- getConditionRunIdMapping(mdQuantExports$proteinIntensitiesLong)
  # intensityConditions <- conditionRunIdMapping[colnames(mdQuantExports$proteinIntensitiesLong)[-1],"Condition"]
  proteinIntensitiesWide <- spreadProteinIntensities(mdQuantExports$proteinIntensitiesLong)

  for (comparison in mdQuantExports$comparisons){

    print(comparison)

    # get the conditions
    comparisonIndex = conditionComparisonMapping$conditionComparison == comparison
    first_condition = conditionComparisonMapping$up.condition[comparisonIndex]
    second_condition = conditionComparisonMapping$down.condition[comparisonIndex]

    print(first_condition)
    print(second_condition)


    # get protein viz data (row data and statistics)
    rowDataStatistics = mdQuantExports$proteinViz[comparisonIndex,]$data[[1]]
    if ("ProteinIds" %in% colnames(rowDataStatistics)){
      print("V2 detected")
      colnames(rowDataStatistics) = c("ProteinGroupId", "GroupLabel", "GroupLabelType", "ProteinId", "GeneName", "Description", "ProteinQValue", "FoldChange", "AdjustedPValue", "PValue", "ConfLow", "ConfHigh")
    } else {
      print("V2 detected")
    }
    rowDataStatistics = handleProteinVizIds(rowDataStatistics)

    # Get protein intensity data (assay data)
    proteinIntensitiesComparison <- filterIntensitiesByProteinViz(rowDataStatistics,
                                                                  proteinIntensitiesWide)
    proteinIntensitiesComparison <- filterProteinIntensityColumnsbyComparison(first_condition,
                                                                              second_condition,
                                                                              conditionRunIdMapping,
                                                                              proteinIntensitiesComparison)

    # Get the column data/experimental design
    columnData <- createColumnData(first_condition,
                                   second_condition,
                                   conditionRunIdMapping)

    comparisonExperimentsList[[comparison]] =
      md2ProteinSummarizedExperiment(rowDataStatistics,
                                     proteinIntensitiesComparison,
                                     columnData,
                                     first_condition,
                                     second_condition,
                                     by = "Protein")
  }

  comparisonExperimentsList
}



#' Get mapping from comparison strings to up and down conditions
#'
#' @param proteinViz a table produced by md quant workflows with DE statistics
#' @return a dataframe containing the up and down conditions for each comparison
#' @export getConditionComparisonMapping
getConditionComparisonMapping <- function(proteinViz){

  stopifnot("up.condition" %in% names(proteinViz))

  conditionComparisonMapping = as.data.frame(cbind(proteinViz$up.condition,
                                                   proteinViz$down.condition,
                                                   proteinViz$conditionComparison),
                                             stringsAsFactors = T)

  colnames(conditionComparisonMapping) = c("up.condition",
                                           "down.condition",
                                           "conditionComparison")

  conditionComparisonMapping
}


#' Long to wide conversion of protein Intensities
#'
#' @param proteinIntensitiesLong a long form protein intensities object produced by the md quant workflow
#' @return a wide dataframe containing all the processed protein intensity information
#' @export spreadProteinIntensities
spreadProteinIntensities <- function(proteinIntensitiesLong){

  #handle TMT imports
  proteinIntensitiesLong = renameTMTIntensityColumn(proteinIntensitiesLong)

  proteinIntensitiesLong = proteinIntensitiesLong[,c("id","log2NInt", "run_id")]
  proteinIntensitiesLong = proteinIntensitiesLong %>% group_by(id, run_id)
  proteinIntensitiesWide = proteinIntensitiesLong %>%
    group_by(id) %>%
    spread(run_id,  log2NInt)

  colnames(proteinIntensitiesWide)[1] = "ProteinGroupId"
  proteinIntensitiesWide
}

#' Returns a dataframe matching Run Id's to Conditions
#' @export getConditionRunIdMapping
getConditionRunIdMapping <- function(proteinIntensitiesLong){
  conditionRunIdMapping <- unique(proteinIntensitiesLong[,c("condition","run_id")])
  colnames(conditionRunIdMapping) <- c("Condition", "IntensityColumn")
  conditionRunIdMapping[] <- lapply(conditionRunIdMapping, as.character)
  rownames(conditionRunIdMapping) <- conditionRunIdMapping$IntensityColumn
  conditionRunIdMapping
}

#' Get's the Intensity Columns for a given Condition
#' @export getIntensityColumnsFromCondition
getIntensityColumnsFromCondition <- function(condition, conditionRunIdMapping){
  IntensityColumns <- conditionRunIdMapping[conditionRunIdMapping$Condition == condition,
                                            "IntensityColumn"]
  IntensityColumns
}



#' filter the protein intensity by columns for a comparison
#' @export filterProteinIntensityColumnsbyComparison
filterProteinIntensityColumnsbyComparison <- function(condition1,
                                                      condition2,
                                                      conditionRunIdMapping,
                                                      assay_data){


  requiredColumns = c("ProteinId", # for protein GroupId
                      getIntensityColumnsFromCondition(condition1, conditionRunIdMapping),
                      getIntensityColumnsFromCondition(condition2, conditionRunIdMapping))
  print(requiredColumns)

  stopifnot(length(requiredColumns) > 2)
  return(assay_data[,requiredColumns])
}

#' Filters protein intensities by the proteins which the limma statistics
#' were calculated upon.
#' @export filterIntensitiesByProteinViz
filterIntensitiesByProteinViz <- function(proteinViz, proteinIntensitiesWide){
  importantColumns = c("ProteinGroupId", "ProteinId", "PValue") # complete.cases for important columns
  completeRows = complete.cases(proteinViz[,importantColumns])
  proteinViz = proteinViz[completeRows,][,c("ProteinGroupId", "ProteinId")]
  proteinIntensitiesWide = merge(proteinViz, proteinIntensitiesWide, by="ProteinGroupId", all.x =T)
  proteinIntensitiesWide[-1]
}

#' Creates the Column Data (experimental design) for a given experimental Comparison
#' @export createColumnData
createColumnData <- function(condition1, condition2, conditionRunIdMapping){

  requiredRows <- as.logical((conditionRunIdMapping$Condition == condition1) +
                               (conditionRunIdMapping$Condition == condition2))
  columnData <- conditionRunIdMapping[which(requiredRows), ]

  columnData$GROUP = (condition2 == columnData$Condition)
  columnData
}

#' This function generates the protein summarized experiment from the processed md artifacts
#' @export md2ProteinSummarizedExperiment
md2ProteinSummarizedExperiment <- function(rowDataStatistics, proteinIntensitiesComparison, columnData, up.condition, down.condition, by = "Protein"){

  # first remove any columns with all NA's
  rowDataStatistics <- rowDataStatistics[complete.cases(rowDataStatistics),]
  proteinIntensitiesComparison <- proteinIntensitiesComparison[complete.cases(proteinIntensitiesComparison),]

  print(paste("Rows in complete protein viz", dim(rowDataStatistics)[1]))
  print(paste("Rows in complete assay data/protein intensity measurement", dim(proteinIntensitiesComparison)[1]))

  nrows <- dim(proteinIntensitiesComparison)[1]
  ncols <- dim(proteinIntensitiesComparison)[2]

  row_data = merge(proteinIntensitiesComparison, rowDataStatistics, by = "ProteinId", all.y = T)
  row_data = row_data[,c("ProteinId","GeneName","PValue","AdjustedPValue","FoldChange")]

  #column name conversions for using EnrichmentBrowser
  colnames(row_data)[3:5] = c("PVAL","ADJ.PVAL","FC")

  proteinIntensitiesComparison$ProteinGroupId <- NULL
  rownames(proteinIntensitiesComparison) <- proteinIntensitiesComparison$ProteinId
  assay_data = as.matrix((proteinIntensitiesComparison[row_data$ProteinId,
                                                       columnData$IntensityColumn]))

  print(dim(columnData))
  print(dim(assay_data))
  rse <- SummarizedExperiment(rowData = row_data,
                              assays= SimpleList(counts=assay_data),
                              colData = columnData)

  metadata(rse) <- list(up.condition = up.condition,
                        down.condition = down.condition)

  if (by == "Protein"){
    ids <- row_data$ProteinId
  } else if (by == "Gene"){
    ids <- toupper(row_data$GeneName)
  }

  rownames(rse) = ids

  # Enforce accurate FC's
  row_data_fc <- rowData(rse)$FC
  calc_fc <- rowMeans(assay(rse)[,rse$GROUP == 0]) - rowMeans(assay(rse)[,rse$GROUP == 1])

  print(sum((row_data_fc - calc_fc)**2)/length(row_data_fc))

  stopifnot((row_data_fc - calc_fc) < 2**(-10))

  rse
}


#' @export convertDiscoveryProteinIntensties
#' @param proteinIntensitiesLong a table coming from the protein int object in the disco workflow
#' @return proteinIntensitiesLong an updated proteinIntensitiesLong that looks like the maxquant workflow output
convertDiscoveryProteinIntensties <- function(proteinIntensitiesLong){
  proteinIntensitiesLong <- proteinIntensitiesLong[,c("run_id",
                                                      "ProteinGroupId",
                                                      "log2NInt_ProteinGroupId",
                                                      "condition")]
  colnames(proteinIntensitiesLong) <- c("run_id",
                                        "id",
                                        "log2NInt",
                                        "condition")
  proteinIntensitiesLong
}

#' @export convertV2ProteinIntensties
#' @param proteinIntensitiesLong a table coming from the protein int object in the disco workflow
#' @return proteinIntensitiesLong an updated proteinIntensitiesLong that looks like the maxquant workflow output
convertV2ProteinIntensties <- function(proteinIntensitiesLong){
  if ('Intensity' %in% colnames(proteinIntensitiesLong)) {
    intensity_col = 'Intensity'
  } else if ('NormalisedIntensity' %in% colnames(proteinIntensitiesLong)) {
    intensity_col = 'NormalisedIntensity'
  }

  proteinIntensitiesLong <- proteinIntensitiesLong[,c("replicate",
                                                      "GroupId",
                                                      intensity_col,
                                                      "condition")]
  colnames(proteinIntensitiesLong) <- c("run_id",
                                        "id",
                                        "log2NInt",
                                        "condition")
  proteinIntensitiesLong$log2NInt = log2(proteinIntensitiesLong$log2NInt)
  proteinIntensitiesLong
}

#' @export handleProteinVizIds
#' @param comparison_viz a table coming from the protein viz object
#' @return comparison_viz the same table with pipe ids handled
handleProteinVizIds <- function(comparison_viz){
  ids <- unlist(lapply(comparison_viz$ProteinId, parsePipeId))
  comparison_viz$ProteinId <- ids
  comparison_viz
}

#' @export parsePipeId
#' @param id one protein id
#' @return parsed protein ID, if Uniprot Only the first protein ID is returned.
parsePipeId <- function(id){

  # assumption that if there is a ; we split by proteins and then check if it's uniprot and take the second element
  if (str_detect(id, ";")){
    split_by_semicolon <- unlist(stringr::str_split(id,";"))
    first_component <- split_by_semicolon[1]

    if(str_detect(first_component, '(sp|tr)\\|(.*)')){
      components <- unlist(stringr::str_split(first_component, '(sp|tr)\\|'))
      id <- components[2]
    }
  }

  if (str_detect(id, 'sp|tr\\|(.*)\\|')){
    components <- unlist(stringr::str_split(id, '\\|'))
    id <- components[2]
  }

  return(id)
}


#' @export renameTMTIntensityColumn
#' @param comparison_viz a table coming from the protein int
#' @return a table with a renamed protein intensity column without "norm"
renameTMTIntensityColumn <- function(proteinIntensitiesLong){

  if ("log2NIntNorm" %in% colnames(proteinIntensitiesLong)){
    proteinIntensitiesLong$log2NInt <- proteinIntensitiesLong$log2NIntNorm
  }

  proteinIntensitiesLong
}
