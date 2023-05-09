#' @export
write_set_infos <- function(list_info_df, outfolder, db_name){
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
  write_csv(list_info_df, file.path(outfolder, paste0(db_name, "_set_info.csv")))
}


#' @import pathwayPCA
#' @export
write_gmt <- function(genesets_proteins, outfolder, db_name){

  list_proteins <- df_items_to_list(genesets_proteins)

  cp_pathwayCollection <- CreatePathwayCollection(
    sets_ls = list_proteins,
    TERMS = names(list_proteins)
  )

  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
  pathwayPCA::write_gmt(
    pathwayCollection = cp_pathwayCollection,
    file = file.path(outfolder, paste0(db_name, ".gmt"))
  )
}


#' Write GMT for GO
#' @export
write_gmt_go <- function(go_sets, outfolder){
  sets_bp <- go_sets$genesets_table$annotation_id[go_sets$genesets_table$annotation_subcategory %in% "biological_process"]
  sets_cp <- go_sets$genesets_table$annotation_id[go_sets$genesets_table$annotation_subcategory %in% "cellular_component"]
  sets_mf <- go_sets$genesets_table$annotation_id[go_sets$genesets_table$annotation_subcategory %in% "molecular_function"]

  proteins_bp <- go_sets$genesets_proteins[go_sets$genesets_proteins$annotation_id %in% sets_bp,]
  proteins_cp <- go_sets$genesets_proteins[go_sets$genesets_proteins$annotation_id %in% sets_cp,]
  proteins_mf <- go_sets$genesets_proteins[go_sets$genesets_proteins$annotation_id %in% sets_mf,]

  EnrichmentAnalysisStepR::write_gmt(genesets_proteins = proteins_bp, outfolder = outfolder, "BiologicalProcess")
  EnrichmentAnalysisStepR::write_gmt(genesets_proteins = proteins_cp, outfolder = outfolder, "CellularComponents")
  EnrichmentAnalysisStepR::write_gmt(genesets_proteins = proteins_mf, outfolder = outfolder, "MolecularFunction")
}

#' Write Set infos for GO
#' @export
write_set_info_go <- function(go_sets, outfolder = "../inst/extdata/"){
  go_sets$genesets_table$go_download_link <- go_sets$versions$download_link

  sets_bp <- go_sets$genesets_table[go_sets$genesets_table$annotation_subcategory %in% "biological_process",]
  sets_bp$annotation_subcategory <- NULL

  sets_cp <- go_sets$genesets_table[go_sets$genesets_table$annotation_subcategory %in% "cellular_component",]
  sets_cp$annotation_subcategory <- NULL

  sets_mf <- go_sets$genesets_table[go_sets$genesets_table$annotation_subcategory %in% "molecular_function",]
  sets_mf$annotation_subcategory <- NULL

  write_set_infos(list_info_df = sets_bp, outfolder = outfolder, db_name = "BiologicalProcess")
  write_set_infos(list_info_df = sets_cp, outfolder = outfolder, db_name = "CellularComponents")
  write_set_infos(list_info_df = sets_mf, outfolder = outfolder, db_name = "MolecularFunction")

}

#' Write Set infos for Reactome
#' @export
write_set_info_reactome <- function(reactome_sets, outfolder = "../inst/extdata/"){
  reactome_sets$genesets_table$go_download_link <- reactome_sets$versions$download_link
  write_set_infos(list_info_df = reactome_sets$genesets_table, outfolder = outfolder, db_name = "Reactome")
}


#' Items table to list format
#' @keywords internal
df_items_to_list <- function(genesets_proteins){
  gmt_list <- genesets_proteins %>%
    split(.$annotation_id) %>%
    lapply(function(x) x$items)
  return(gmt_list)
}
