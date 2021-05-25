# post process output, format nicely

#' @param results a results table producing output
#' @return The input table with an added column ADJ.PVAL with BH adjusted pvalues.
#' @examples
#'
#' @export enrichment_adjust_p_values
enrichment_adjust_p_values <- function(results){
  for (condition in names(results)){
    for (library in names(results[[condition]])){
      p = results[[condition]][[library]]$PVAL
      results[[condition]][[library]]$ADJ.PVAL = p
    }
  }

  results
}

#' @param results a results table producing output
#' @param gmt_folder a folder which should contain annotation csv's with details for each gmt used in enrichment
#' @return The input table with a added columns with gene set information.
#' @examples
#'
#' @export enrichment_annotate_results
enrichment_annotate_results <- function(results, gmt_folder){
  for (condition in names(results)){
    for (library in names(results[[condition]])){
      enr = results[[condition]][[library]]
      info_file = paste(gsub(".gmt","",library),"_set_info.csv", sep ="")
      info = read.csv(file.path(gmt_folder,info_file), stringsAsFactors = F)
      colnames(info) = c("GENE.SET", "NAME", "DESCRIPTION", "ITEMS", "SIZE")
      if ("ITEMS" %in% colnames(info)){
        info$ITEMS <- NULL
      }
      
      print(colnames(enr))
      print(colnames(info))
      
      enr = as.data.frame(merge(enr, info))
      colnames(enr) <- tolower(colnames(enr))
      results[[condition]][[library]] = enr
    }
  }

  results
}


#' @param results a results table producing output
#' @return Writes json objects with enrichment results into the output folder
#' @examples
#'
#' @export enrich_write_output_tables
enrich_write_output_tables <- function(results, output_folder){
  for (comparison in names(results)){
    first_condition = strsplit(comparison,"\\s+")[[1]][1]
    second_condition = strsplit(comparison,"\\s+")[[1]][3]
    for (library in names(results[[comparison]])){
      file_name = str_c(str_c(first_condition, second_condition,sep = "."),
                        library,sep = ".")

      enr = list(up.condition = unbox(first_condition),
                 down.condition = unbox(second_condition),
                 database = unbox(library),
                 version = unbox("0.0.0"),
                 gene.set.statistics = as.data.frame(results[[comparison]][[library]]))
      
      print(str_c(file_name,".json"))
      write_json(enr,file.path(output_folder,str_c(file_name,".json")))
    }
  }
}

# utility for fixing strings in columns containing arrays from read_csv
.string_to_array <- function(bad_string){
  bad_string <- gsub("'","",bad_string)
  bad_string <- gsub("\\[","",bad_string)
  bad_string <- gsub("\\]","",bad_string)
  vec <- str_split(bad_string, ", ")
  vec
}