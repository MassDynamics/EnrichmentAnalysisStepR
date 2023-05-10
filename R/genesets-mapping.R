
#' Map sets to human R files
#'
#' @param set Name of the annotation set
#' @export
map_sets_to_file_human <- function(set){
  file <- switch(
    set,
    "BiologicalProcess" = "GO_human_2023-05-08.rda",
    "CellularComponents" = "GO_human_2023-05-08.rda",
    "MolecularFunction" = "GO_human_2023-05-08.rda",
     "Reactome"= "msigdb_v7.5.rda",
     "MsigDB-C2.CGP"= "msigdb_v7.52.rda",
     "MsigDB-C2.CP"= "msigdb_v7.52.rda",
     "MsigDB-C3"= "msigdb_v7.52.rda",
     "MsigDB-C4"= "msigdb_v7.52.rda",
     "MsigDB-C5.HPO"= "msigdb_v7.52.rda",
     "MsigDB-C6"= "msigdb_v7.52.rda",
     "MsigDB-C7"= "msigdb_v7.52.rda",
     "MsigDB-C8"= "msigdb_v7.52.rda",
     "MsigDB-h.all." = "msigdb_v7.52.rda"
  )
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
    "BiologicalProcess" = "GO_mouse_2023-05-08.rda",
    "CellularComponents" = "GO_mouse_2023-05-08.rda",
    "MolecularFunction" = "GO_mouse_2023-05-08.rda",
    "Reactome"= "Reactome_2023-05-08.rda",
  )
}

#' Map sets to yeast R files
#'
#' @param set Name of the annotation set
#' @export
map_sets_to_file_yeast <- function(set){
  file <- switch(
    set,
    "BiologicalProcess" = "GO_yeast_2023-05-08.rda",
    "CellularComponents" = "GO_yeast_2023-05-08.rda",
    "MolecularFunction" = "GO_yeast_2023-05-08.rda",
    "Reactome"= "Reactome_2023-05-08.rda",
  )
}
