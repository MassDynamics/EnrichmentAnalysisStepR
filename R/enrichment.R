#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#'
#' @param comparison A summarized experiment object containing only two experimental groups
#' @param gmt_folder a folder containing gene set libraries which need to be referenced
#' @param method the method enrichment browser should use for set based enrichment
#' @return enrichment An enrichment table provided by EnrichmentBrowserPackage
#' @examples
#' sdgdfsgsf
#' @export perform_comparison_enrichment
perform_comparison_enrichment <- function(comparison, gmt_folder, method){

  enrichment = list()

  for (library in list.files(gmt_folder, pattern = "\\.gmt$")){
    gmt.file <- file.path(gmt_folder, library)

    gs <- getGenesets(gmt.file)

    print(gmt.file)
    print(paste("Sets in Library: ", length(gs)))
    print(paste("Genes in Library: ", length(unique(unlist(gs)))))
    print(paste("Genes in Library and dataset: ", length(intersect(rownames(comparison),unlist(gs)))))

    # check for common ids
    if (length(intersect(rownames(comparison),unlist(gs)))>100){

      # do enrichment
      sbea.res <- sbea(method = method, se = comparison, gs = gmt.file, perm = 100, alpha = 1)
      results = gsRanking(sbea.res)

      # get sets info (size/observed and merge)
      #sets_info <- add_sets_info(gs, comparison)
      #results <- merge(results,sets_info)

      # write results
      lib = str_split(library,".gmt")[[1]][1]
      label = str_c(lib,"_enrichment.csv")

      enrichment[[lib]] = results
    } else {
      print("Skipping enrichment, not enough genes in common.")
      print("Heres some example ids")
      print(rownames(comparison)[1:10])
      print(unlist(gs)[1:10])
    }

  }

  enrichment
}

#' @param upload_folder a folder containing protein_viz and protein_counts and intensity files
#' @param gmt_folder a folder containing gene set libraries which need to be referenced
#' @param output_folder a folder where the output files are deposited
#' @param by either "Gene" or "Protein" which will determine whether genes or proteins are used in the gene set libraries.
#' @param method the method enrichment browser should use for set based enrichment
#' @return Writes enrichment outputs to the output folder.
#' @examples
#' enrichment_workflow_step(upload_folder,
#'  gmt_folder = gmt_folder,
#'  output_folder = output_folder,
#'  by = by)
#' @export enrichment_workflow_step
enrichment_workflow_step <- function(upload_folder,
                                     gmt_folder ="./ckg_gmts",
                                     output_folder = "./gsea_results",
                                     by = "Protein",
                                     method = "gsea"){


  # Get Experiment Data
  protein_viz_path = file.path(upload_folder,"protein_viz.json")
  stopifnot(file.exists(protein_viz_path))
  protein_viz = read_json(protein_viz_path, simplifyVector = TRUE)
  stopifnot(dim(protein_viz)[1]>0)

  protein_count_and_intensities = read_json(file.path(upload_folder,"protein_counts_and_intensity.json"), simplifyVector = TRUE)
  protein_count_and_intensities = process_protein_intensities(protein_count_and_intensities, by = by)
  stopifnot(dim(protein_count_and_intensities)[1]>0)

  # for comparison in comparisons...

  results = list()
  for (row in (1:nrow(protein_viz))){

    comparison = protein_viz[row,"conditionComparison"]
    print(comparison)

    comparison_viz = protein_viz[row, "data"]
    tmp_protein_int = filter_protein_ints(comparison, protein_count_and_intensities)

    cls_vec = get_cls_vec(comparison, tmp_protein_int)
    se_comparison <- protein_viz_int_2_de_exp(comparison_viz, tmp_protein_int, cls_vec, by = by)

    results[[comparison]] = perform_comparison_enrichment(se_comparison, gmt_folder, method)

  }

  # compile output tables by merging enrichment results with annotation data
  results <- enrichment_adjust_p_values(results)
  results <- enrichment_annotate_results(results, gmt_folder)
  enrich_write_output_tables(results, output_folder)


  # done!
  results
}


add_sets_info <-function(gs, se){

  set_membership = stack(gs)
  set_size = as.data.frame(table(set_membership$ind))

  sets_info = set_membership %>% group_by(ind) %>%
    summarize(genes = list(sort(unique(values))))

  sets_info$set_size = unlist(lapply(sets_info$genes, length))
  sets_info$set_intersection = unlist(lapply(sets_info$genes, function (x) length(intersect(x,rownames(se)))))
  colnames(sets_info)
  colnames(sets_info) <- c("GENE.SET", "ITEMS", "SIZE", "OBSERVED")

  sets_info
}
