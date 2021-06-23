#' Automates Enrichment Analyses on Summarized Objects for a library of gene sets
#'
#' @param comparison A summarized experiment object containing only two experimental groups
#' @param gmt_folder a folder containing gene set libraries which need to be referenced
#' @param method the method enrichment browser should use for set based enrichment
#' @return enrichment An enrichment table provided by EnrichmentBrowserPackage
#' @examples
#' sdgdfsgsf
#' @export perform_comparison_enrichment
perform_comparison_enrichment <- function(summarized_experiment, gmt_folder, method, perm = perm){
  
  enrichment = list()
  
  for (library in list.files(gmt_folder, pattern = "\\.gmt$")){
    gmt.file <- file.path(gmt_folder, library)
    
    gs <- getGenesets(gmt.file)
    
    print(gmt.file)
    print(paste("Sets in Library: ", length(gs)))
    print(paste("Genes in Library: ", length(unique(unlist(gs)))))
    print(paste("Genes in Library and dataset: ", length(intersect(rownames(summarized_experiment),unlist(gs)))))
    
    # check for common ids
    if (length(intersect(rownames(summarized_experiment),unlist(gs)))>100){
      
      # do enrichment      
      if (method == "camera"){
        results <- .camera(summarized_experiment, gs)
      } else {
        sbea.res <- sbea(method = method, se = summarized_experiment, 
                         gs = gmt.file, perm = perm, alpha = 1)
        results = as.data.frame(gsRanking(sbea.res))
        colnames(results) <- tolower(colnames(results))
      }
      
      print(colnames(results))
      # pval and gene.set must be present!
      for (col in c("gene.set", "pval")){
        stopifnot(col %in% colnames(results))
      }
      
      # write results
      lib = str_split(library,".gmt")[[1]][1]
      label = str_c(lib,"_enrichment.csv")
      
      results <- add_gene_sets_to_result(results, gs)
      results <- add_AFC_to_results(results, summarized_experiment)
      
      # add adjusted p value
      results$adj.pval = p.adjust(results$pval, method = "BH") 
      
      print(colnames(results))
      stopifnot("afc" %in% colnames(results))
      stopifnot("items" %in% colnames(results))
      
      enrichment[[lib]] = results
      
    } else {
      print("Skipping enrichment, not enough genes in common.")
      print("Heres some example ids")
      print(rownames(summarized_experiment)[1:10])
      print(unlist(gs)[1:10])
    }
    
  }
  
  enrichment
}




.camera <- function(summarized_experiment, gene_sets){
  #Part of this code was derived from the DiscoveryBrowser package. 
  
  
  print("Running CAMERA using native code")
  
  
  # Prepare CAMERA Inputs
  expression_data <- assay(summarized_experiment)
  grp <- colData(summarized_experiment)$GROUP
  group <- factor(grp)
  f <- "~" 
  f <- formula(paste0(f, "group"))
  design <- model.matrix(f)
  
  # filter for genes in the experiment
  igenes <- intersect(rownames(expression_data), unique(unlist(gene_sets)))
  if(!length(igenes)) stop("Expression dataset (se)", " and ", 
                           "gene sets (gs) have no gene IDs in common")
  expression_data <- expression_data[igenes,]
  filtered_gene_sets <- lapply(gene_sets, function(s) s[s %in% igenes])
  
  
  # Call CAMERA
  gs_index <- limma::ids2indices(filtered_gene_sets, rownames(expression_data))
  camera_result <- limma::camera(expression_data, gs_index, design, sort=FALSE)
  
  # format results
  camera_result$Pathway <- rownames(camera_result)
  
  camera_result <- camera_result[, c("Pathway","PValue", "FDR","NGenes")]
  colnames(camera_result) <- c("gene.set", "pval", "adj.pval", "observed")
  
  camera_result
}

add_gene_sets_to_result <- function(result_table, gs){
  result_table$items <- gs[match(result_table$gene.set, names(gs))]
  result_table
}

# utility for converting list of char vectors to table
#' @export gs2df
gs2df <- function(gs){
  df <- as.data.frame(cbind(names(gs),gs))
  colnames(df) <- c("gene.set","items")
  df$gene.set <- sapply(df$gene.set, toString)
  df
}

# get's average fold change using a gene set and the SE artifact.
#' @export get_AFC
get_AFC <- function(summarized_experiment, gene_set){
  mean(rowData(summarized_experiment)[rownames(summarized_experiment) %in% gene_set,"FC"])
}

#' @export add_AFC_to_results
add_AFC_to_results <- function(results, summarized_experiment){
  results$afc <- unlist(lapply(results$items, function(x){get_AFC(summarized_experiment, x)}))
  results
}
