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
                                     method = "gsea",
                                     perm = 250){
  
  
  # Get Experiment Data
  protein_viz_path = file.path(upload_folder,"protein_viz.json")
  stopifnot(file.exists(protein_viz_path))
  protein_viz = read_json(protein_viz_path, simplifyVector = TRUE)
  stopifnot(dim(protein_viz)[1]>0)
  
  
  protein_int_path = file.path(upload_folder,"proteinGroups_quant_intensities.txt")
  if (!file.exists(protein_int_path)){
    protein_int_path = file.path(upload_folder,"Protein_Intensity.txt")
    stopifnot(file.exists(protein_int_path))
    protein_int = read.csv(protein_int_path, sep ="\t")
    protein_int = disco_protein_int_to_byo_int(protein_int)
  } else {
    stopifnot(file.exists(protein_int_path))
    protein_int = read.csv(protein_int_path, sep ="\t")
  }
  
  protein_int <- get_protein_quant_intensities(protein_int)
  
  # for comparison in comparisons...
  results = list()
  for (row in (1:nrow(protein_viz))){
    
    # Get the Comparison
    comparison = protein_viz[row,"conditionComparison"]
    print(comparison)
    
    # get protein viz data
    comparison_viz = as.data.frame(protein_viz[row, "data"])
    comparison_viz = handle_protein_viz_ids(comparison_viz)
    
    #Get the rest of the data
    assay_data <- filter_int_by_viz(comparison_viz, protein_int)
    assay_data <- filter_prot_int_cols_by_comp(comparison, assay_data)
    cls_vec <- get_cls_vec(comparison, assay_data)
    se_comparison <- md_to_summarized_experiment(comparison_viz, assay_data, cls_vec, by = by)
    
    results[[comparison]] = perform_comparison_enrichment(se_comparison, gmt_folder, method, perm)
    
  }
  
  # compile output tables by merging enrichment results with annotation data
  results <- enrichment_adjust_p_values(results)
  results <- enrichment_annotate_results(results, gmt_folder)
  enrich_write_output_tables(results, output_folder)
  
  
  # done!
  results
}
