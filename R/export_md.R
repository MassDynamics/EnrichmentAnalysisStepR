
# MD Specific JSON output writing and annotation integration

#' @param listEnrichmentResults a listEnrichmentResults table producing output
#' @param gmt_folder a folder which should contain annotation csv's with details for each gmt used in enrichment
#' @return The input table with a added columns with gene set information.
#' @examples
#'
#' @export annotateEnrichmentResults
annotateEnrichmentResults <- function(listEnrichmentResults, gmt_folder){

  for (comparison in names(listEnrichmentResults)){
    annotateComparison(comparison, gmt_folder)
  }

  listEnrichmentResults
}

annotateComparison <- function(comparison, gmt_folder){
  for (library in names(comparison$libraryStatistics)){
    enr = comparison$libraryStatistics[[library]]

    info_file = paste(gsub(".gmt","",library),"_set_info.csv", sep ="")
    if (file.exists(file.path(gmt_folder,info_file))){
      print(paste("Using set info file: ", info_file))
      info = read.csv(file.path(gmt_folder,info_file), stringsAsFactors = F)
      colnames(info) = c("gene.set", "name", "description", "items", "size", "linkout")
      if ("items" %in% colnames(info)){
        info$items <- NULL
      }

      enr = as.data.frame(merge(enr, info), by = "gene.set")

    } else {
      print(paste("Could not file info file: ", file.path(gmt_folder,info_file)))
      enr = as.data.frame(enr)
    }
    colnames(enr) <- tolower(colnames(enr))
    comparison$libraryStatistics[[library]] = enr
  }
  return(comparison)
}


#' @param listEnrichmentResults a listEnrichmentResults table producing output
#' @return Writes json objects with enrichment listEnrichmentResults into the output folder
#' @examples
#'
#' @export writeEnrichmentJSON
writeEnrichmentJSON <- function(listEnrichmentResults, outputFolder){

  dir.create(outputFolder, showWarnings = FALSE)

  for (comparison in names(listEnrichmentResults)){
    print(paste0("Write results for comparison: ", comparison))
    if(is.null(names(listEnrichmentResults[[comparison]]$libraryStatistics))){
      print(paste0("No enrichment available for comparison: ", comparison))
    }
    for (library in names(listEnrichmentResults[[comparison]]$libraryStatistics)){
      print(paste0("Inside for:", comparison))
      statistics = listEnrichmentResults[[comparison]]$libraryStatistics[[library]]
      up.condition = listEnrichmentResults[[comparison]]$up.condition
      down.condition = listEnrichmentResults[[comparison]]$down.condition

      enr = list(up.condition = unbox(up.condition),
                 down.condition = unbox(down.condition),
                 database = unbox(library),
                 version = unbox("0.0.0"),
                 gene.set.statistics = as.data.frame(statistics))

      # remove special characters from file name
      file_name = str_c(comparison,library,sep = ".")
      file_name = gsub("[[:punct:]]", "_", file_name)
      file_name = gsub(" ", "_", file_name)
      print(str_c(file_name,".json"))
      write_json(enr,file.path(outputFolder,str_c(file_name,".json")), digits = NA)
    }
  }
}
