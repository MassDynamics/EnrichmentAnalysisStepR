# md experiment parsing


#' @param protein_viz a massdynamics protein visualization json that has been read into a dataframe
#' @param protein_ints a processed protein intensity dataframe produced by process protein intensities
#' @param cls_vec a class vector indicating the classes of columns used to
#' @return a summarized experiment object with row data for differential expression and column intensity data
#' @examples
#'
#' @export md_to_summarized_experiment
md_to_summarized_experiment <- function(protein_viz, protein_ints, cls_vec, by = "Protein"){
  
  # remove data where we didn't impute results (missing complete condition)
  
  # first remove any columns with all NA's
  experiment <- protein_ints[,colSums(is.na(protein_ints))<nrow(protein_ints)]
  experiment <- protein_ints[complete.cases(protein_ints),]
  
  cls = colnames(protein_ints)
  cls = cls[2:length(cls)]
  
  
  nrows <- dim(experiment)[1]
  ncols <- dim(experiment)[2]
  
  md_colData <- DataFrame(Treatment=cls)
  
  row_data = merge(experiment, protein_viz)#, by = "ProteinId", all.y = T)
  row_data = row_data[,c("ProteinId","GeneName","PValue","AdjustedPValue","FoldChange")]
  colnames(row_data)[3:5] = c("PVAL","ADJ.PVAL","FC")
  
  experiment$ProteinGroupId <- NULL
  assay_data = as.matrix((experiment[,cls]))
  
  
  rse <- SummarizedExperiment(rowData = row_data,
                              assays=SimpleList(counts=assay_data),
                              colData = md_colData)
  
  if (by == "Protein"){
    ids <- row_data$ProteinId
  } else if (by == "Gene"){
    ids <- toupper(row_data$GeneName)
  }
  
  rownames(rse) = ids
  rse$GROUP <- as.numeric(cls_vec)
  rse
}


#' @param comparison a string with a " - " seperator delineating a comparison on the left and right
#' @param protein_counts_and_intensities a protein intensity dataframe produced by read the mass dynamics file
#' @return a processed protein intensity dataframe with columns only for that comparison
#' @examples
#'
#' @export filter_prot_int_cols_by_comp
filter_prot_int_cols_by_comp <- function(comparison, protein_count_and_intensities){
  
  cols = unlist(gsub("_[0-9]*$", "", colnames(protein_count_and_intensities)))
  
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

#' @export filter_int_by_viz
filter_int_by_viz <- function(protein_viz, protein_int){
  important_columns = c("ProteinGroupId", "ProteinId", "PValue") # complete.cases for important columns
  complete_rows = complete.cases(protein_viz[,important_columns])
  protein_viz_comparison = protein_viz[complete_rows,][,c("ProteinGroupId", "ProteinId")]
  assay_data = merge(protein_viz_comparison, protein_int, by="ProteinGroupId", all.x =T)
  assay_data[-1]
}


#' @export get_cls_vec
get_cls_vec <- function(comparison, tmp_protein_int){

  second_condition = strsplit(comparison,"\\s+")[[1]][3]
  cols = unlist(gsub("_[0-9]*$", "", colnames(tmp_protein_int)))
  cls_vec <- grepl(second_condition, cols)
  cls_vec <- cls_vec[-1]

  #handle bigger table
  if ("GeneName" %in% colnames(tmp_protein_int)){
    cls_vec <- cls_vec[-1]
  }
  
  if ("ProteinGroupId" %in% colnames(tmp_protein_int)){
    cls_vec <- cls_vec[-1]
  }

  print(length(cls_vec))
  cls_vec
}


#' @export get_protein_quant_intensities
get_protein_quant_intensities <- function(protein_int){
  
  protein_int = protein_int[,c("id", "condition","log2NInt", "run_id")]
  protein_int = protein_int %>% group_by(id, condition, run_id)
  
  protein_int = protein_int %>%
    group_by(id) %>%
    unite("ID_type", condition, run_id,remove = TRUE) %>%
    spread(ID_type,  log2NInt)
  
  colnames(protein_int)[1] = "ProteinGroupId"
  
  protein_int
}

#' @export handle_protein_viz_ids
#' @param comparison_viz a table coming from the protein viz object
#' @return comparison_viz the same table with pipe ids handled
handle_protein_viz_ids <- function(comparison_viz){
  ids <- unlist(lapply(comparison_viz$ProteinId, parse_pipe_id))
  comparison_viz$ProteinId <- ids
  comparison_viz
}

parse_pipe_id <- function(id){
  
  if (length((grep("\\|", id))>0)){
    components <- unlist(stringr::str_split(id,"\\|"))
    id <- components[2]
  }
  
  return(id)
}
