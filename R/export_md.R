# post process output, format nicely

#' @param results a results table producing output
#' @param gmt_folder a folder which should contain annotation csv's with details for each gmt used in enrichment
#' @return The input table with a added columns with gene set information.
#' @examples
#'
#' @export enrichment_annotate_results
enrichment_annotate_results <- function(results, gmt_folder){
  for (condition in names(results)){
    for (library in names(results[[condition]]$enrichmentResults)){
      enr = results[[condition]]$enrichmentResults[[library]]
      
      
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
      results[[condition]]$enrichmentResults[[library]] = enr
    }
  }
  
  results
}


#' @param results a results table producing output
#' @return Writes json objects with enrichment results into the output folder
#' @examples
#'
#' @export enrich_write_output_tables
enrich_write_output_tables <- function(results, output_folder, conditions){
  
  dir.create(output_folder, showWarnings = FALSE)
  
  for (comparison in names(results)){
    
    first_condition <- results[[comparison]]$up.condition
    second_condition <- results[[comparison]]$down.condition
    
    for (library in names(results[[comparison]]$enrichmentResults)){
      file_name = str_c(str_c(first_condition, second_condition,sep = "."),
                        library,sep = ".")
      
      enr = list(up.condition = unbox(first_condition),
                 down.condition = unbox(second_condition),
                 database = unbox(library),
                 version = unbox("0.0.0"),
                 gene.set.statistics = as.data.frame(results[[comparison]]$enrichmentResults[[library]]))
      
      # remove special characters from file name 
      file_name = gsub("[[:punct:]]", "", file_name) 
      print(str_c(file_name,".json"))
      write_json(enr,file.path(output_folder,str_c(file_name,".json")), digits = NA)
    }
  }
}