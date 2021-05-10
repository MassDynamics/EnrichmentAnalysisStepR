# md experiment parsing


#' @param protein_viz a massdynamics protein visualization json that has been read into a dataframe
#' @param protein_ints a processed protein intensity dataframe produced by process protein intensities
#' @param cls_vec a class vector indicating the classes of columns used to
#' @return a summarized experiment object with row data for differential expression and column intensity data
#' @examples
#'
#' @export protein_viz_int_2_de_exp
protein_viz_int_2_de_exp <-function(protein_viz, protein_ints, cls_vec, by = "Protein"){

  # remove data where we didn't impute results (missing complete condition)
  experiment <- protein_ints[complete.cases(protein_ints),]
  cls = gsub("_[0-9]*","",colnames(protein_ints))

  if (by == "Protein"){
    ids <- experiment$ProteinId
    cls = cls[2:length(cls)]
  } else if (by == "Gene"){
    ids <- toupper(experiment$GeneName)
    cls = cls[3:length(cls)]
  }


  nrows <- dim(experiment)[1]
  ncols <- dim(experiment)[2]

  md_colData <- DataFrame(Treatment=cls)

  row_data = merge(experiment, protein_viz, by = "ProteinId")
  row_data = row_data[,c(8,9,14,15,16)]
  colnames(row_data)[3:5] = c("PVAL","ADJ.PVAL","FC")


  if (by == "Protein"){
    assay_data = as.matrix((experiment[,2:ncols]))
  } else if (by == "Gene"){
    assay_data = as.matrix((experiment[,3:ncols]))
  }

  rse <- SummarizedExperiment(rowData = row_data,
                              assays=SimpleList(counts=assay_data),
                              colData = md_colData)

  rownames(rse) = ids
  rse$GROUP <- as.numeric(cls_vec)
  rse
}


#' @param comparison a string with a " - " seperator delineating a comparison on the left and right
#' @param protein_counts_and_intensities a protein intensity dataframe produced by read the mass dynamics file
#' @return a processed protein intensity dataframe with columns only for that comparison
#' @examples
#'
#' @export filter_protein_ints
filter_protein_ints <- function(comparison, protein_count_and_intensities){

  cols = unlist(lapply(strsplit(colnames(protein_count_and_intensities),"_"), "[",1))

  print(comparison)
  first_condition = strsplit(comparison,"\\s+")[[1]][1]
  second_condition = strsplit(comparison,"\\s+")[[1]][3]

  required_columns = grepl(first_condition, cols) + grepl(second_condition, cols)
  required_columns[1] <- 1

  if ("GeneName" %in% colnames(protein_count_and_intensities)){
    print("Allowing gene name")
    required_columns[2] <- 1
  }

  tmp_protein_int = protein_count_and_intensities[, which(1==(required_columns))]

  tmp_protein_int
}


#' @export get_cls_vec
get_cls_vec <- function(comparison, tmp_protein_int){

  second_condition = strsplit(comparison,"\\s+")[[1]][3]
  cols <- unlist(lapply(strsplit(colnames(tmp_protein_int),"_"), "[",1))
  cls_vec <- grepl(second_condition, cols)
  cls_vec <- cls_vec[-1]

  #handle bigger table
  if ("GeneName" %in% colnames(tmp_protein_int)){
    cls_vec <- cls_vec[-1]
  }

  print(length(cls_vec))
  cls_vec
}

#' @export process_protein_intensities
process_protein_intensities <- function(protein_count_and_intensities, by = "Protein"){
  if (by == "Protein"){
    protein_count_and_intensities = process_protein_intensities_by_protein(protein_count_and_intensities)
  }
  else if(by == "Gene"){
    protein_count_and_intensities = process_protein_intensities_by_gene(protein_count_and_intensities)
  }
}

#' @export process_protein_intensities_by_protein
process_protein_intensities_by_protein <- function(protein_count_and_intensities){

  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(conditions)) #first level unwrap
  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(intensityValues)) #second level unwrap
  protein_count_and_intensities = protein_count_and_intensities[,c(2,8,14)]
  protein_count_and_intensities = protein_count_and_intensities %>% group_by(ProteinId,name)  %>% mutate(replicate = row_number())
  protein_count_and_intensities = protein_count_and_intensities %>%
    group_by(ProteinId) %>%
    unite("ID_type", name, replicate,remove = TRUE) %>%
    spread(ID_type,  log2NInt_ProteinGroupId)

  protein_count_and_intensities
}

#' @export process_protein_intensities_by_gene
process_protein_intensities_by_gene <- function(protein_count_and_intensities){

  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(conditions)) #first level unwrap
  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(intensityValues)) #second level unwrap
  print(colnames(protein_count_and_intensities))
  protein_count_and_intensities = protein_count_and_intensities[,c(2,3,8,14)]
  protein_count_and_intensities = protein_count_and_intensities %>% group_by(GeneName, ProteinId,name)  %>% mutate(replicate = row_number())
  protein_count_and_intensities = protein_count_and_intensities %>%
    group_by(GeneName, ProteinId) %>%
    unite("ID_type", name, replicate,remove = TRUE) %>%
    spread(ID_type,  log2NInt_ProteinGroupId)

  protein_count_and_intensities
}
