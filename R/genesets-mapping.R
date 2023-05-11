
#' Map sets to human R files
#'
#' @param set Name of the annotation set
#' @export
map_sets_to_file_human <- function(set){
  file <- switch(
    set,
    "BiologicalProcess" = "GO_MolecularFunction.rda",
    "CellularComponents" = "GO_CellularComponents.rda",
    "MolecularFunction" = "GO_BiologicalProcess.rda",
     "Reactome"= "msigdb_v7.5.rda",
     "MsigDB-C2.CGP"= "msigdb_v7.5.rda",
     "MsigDB-C2.CP"= "msigdb_v7.5.rda",
     "MsigDB-C3"= "msigdb_v7.5.rda",
     "MsigDB-C4"= "msigdb_v7.5.rda",
     "MsigDB-C5.HPO"= "msigdb_v7.5.rda",
     "MsigDB-C6"= "msigdb_v7.5.rda",
     "MsigDB-C7"= "msigdb_v7.5.rda",
     "MsigDB-C8"= "msigdb_v7.5.rda",
     "MsigDB-h.all." = "msigdb_v7.5.rda"
  )
  return(file)
}

#' Map sets to mouse R files
#'
#' @param set Name of the annotation set
#'
#' @return
#' @export
map_sets_to_file_mouse <- function(set){
  file <- switch(
    set,
    "BiologicalProcess" = "GO_mouse.rda",
    "CellularComponents" = "GO_mouse.rda",
    "MolecularFunction" = "GO_mouse.rda",
    "Reactome"= "Reactome.rda",
  )
  return(file)
}

#' Map sets to yeast R files
#'
#' @param set Name of the annotation set
#' @export
map_sets_to_file_yeast <- function(set){
  file <- switch(
    set,
    "BiologicalProcess" = "GO_yeast.rda",
    "CellularComponents" = "GO_yeast.rda",
    "MolecularFunction" = "GO_yeast.rda",
    "Reactome"= "Reactome.rda",
  )
  return(file)
}

#' Shared column names between annotation sets
#' @keywords internal
getSharedColNamesAnnotations <- function(){
  return(c("annotation_name", "annotation_id", "linkout","annotation_description","size"))
}

#' Map GO names
#' @keywords internal
mapGOCategoryNames <- function(lowerCaseName){
  capitalCaseName <- switch(lowerCaseName,
                            "biological_process" = "BiologicalProcess" ,
                            "cellular_component" = "CellularComponents",
                            "molecular_function" = "MolecularFunction")
  return(capitalCaseName)
}
