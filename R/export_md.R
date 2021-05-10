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
      info = read.csv(file.path(gmt_folder,info_file))

      # remove linkout for now
      if ("linkout" %in% colnames(info)){
        colnames(info) <- c("GENE.SET", "LINKOUT","NAME", "DESCRIPTION", "ITEMS", "SIZE")
      } else {
        colnames(info) <- c("GENE.SET", "NAME", "DESCRIPTION", "ITEMS", "SIZE")
      }


      enr = as.data.frame(merge(enr, info))
      results[[condition]][[library]] = enr
    }
  }

  results
}


#' @param results a results table producing output
#' @return Write json objects with enrichment results into the output folder
#' @examples
#'
#' @export enrich_write_output_tables
enrich_write_output_tables <- function(results, output_folder){
  for (condition in names(results)){
    out = file.path(output_folder, condition)
    dir.create(out, recursive = TRUE)
    for (library in names(results[[condition]])){
      enr = as.data.frame(results[[condition]][[library]])
      write_json(enr,file.path(out,str_c(library,".json")))
    }
  }
}
