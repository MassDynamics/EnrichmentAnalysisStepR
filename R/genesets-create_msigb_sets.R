# Some code adapted from: https://github.com/igordot/msigdbr/blob/main/data-raw/msigdbr-prepare.R

#' generate_msig_sets
#'
#' @param xml_filepath file path to msigDB xml file
#' @param uniprot_table_path File path to table downloaded from uniprot
#' @param msig_version msig version
#'
#' @import glue
#' @import dplyr
#' @export generate_msig_sets
generate_msig_sets <- function(xml_filepath,
                               msig_version,
                               uniprot_table_path,
                               max_prop_unmapped = 0.8){

  log4r::info(EnrichmentLogger(), glue("MsigDB version: {msig_version}"))
  log4r::info(EnrichmentLogger(), glue("Uniprot version: {uniprot_table_path}"))
  versions <- data.frame(annotation_name = c("MsigDB","Uniprot"),
                         download_link = c(NA,NA),
                         version = c(msig_version, uniprot_table_path))

  log4r::info(EnrichmentLogger(), "Read and parse MsigDB xml file.")
  mdb_tbl <- msigdb_xml_to_table(xml_filepath)

  log4r::info(EnrichmentLogger(), "Initialise geneset table")
  genesets_table <- create_gene_set_table(mdb_tbl)

  log4r::info(EnrichmentLogger(), "Map gene names to uniprot ids")
  genesets_items <- map_items_to_uniprot(mdb_tbl, uniprot_table_path)
  genesets_items_unmapped <- compute_unmapped_prop(genesets_items)

  log4r::info(EnrichmentLogger(), "Add proportion of unmapped genes to gene sets table.")
  genesets_table <- genesets_table |> dplyr::left_join(genesets_items_unmapped)
  remove_set_ids <- define_genesets_to_remove(genesets_table, max_prop_unmapped)
  log4r::info(EnrichmentLogger(), glue("N = {length(remove_set_ids)} gene sets has > 80% of the genes that could not be mapped. These sets are excluded from the final protein sets."))

  genesets_items_final <- genesets_items |>
    dplyr::filter(!(gs_id %in% remove_set_ids), !is.na(UNIPROTID))
  genesets_items_final <- genesets_items_final |> dplyr::select(gs_id, UNIPROTID) |> dplyr::rename(items = UNIPROTID)

  genesets_table_final <- genesets_table |> dplyr::filter(!(gs_id %in% remove_set_ids))

  log4r::info(EnrichmentLogger(), "Rename output columns")
  genesets_items_final <- rename_genesets_items_columns(genesets_items_final)
  genesets_table_final <- rename_genesets_table_columns(genesets_table_final)
  genes_protein_mapping <- rename_genes_protein_mapping_columns(genesets_items)

  check_genesets_left(genesets_items_final=genesets_items_final,
                      genesets_table_final=genesets_table_final,
                      msig_table=mdb_tbl,
                      max_prop_unmapped=max_prop_unmapped,
                      sets_removed = remove_set_ids)

  log4r::info(EnrichmentLogger(), "Complete.")

  msig_sets <- list(genesets_table = genesets_table_final,
       genesets_proteins = genesets_items_final,
       versions = versions,
       genes_protein_mapping = genesets_items)


  msigSets <- new("AnnotationSetsMsigDB", genesets_table=msig_sets$genesets_table,
                          genesets_proteins=msig_sets$genesets_proteins,
                          versions=msig_sets$versions,
                  genes_protein_mapping=msig_sets$genes_protein_mapping)
}


#' compute_unmapped_prop
#'
#' @param geneset_genes table with gene sets and proteins/genes mapped to them
#' @param genesets_table table with gene sets and their infos
#' @import dplyr
check_genesets_left <- function(genesets_items_final, genesets_table_final, msig_table, max_prop_unmapped, sets_removed){
  # subset of pairs of gs_name/gs_id should be the same in the initial and in the final

  # check that the filtering for prop_unmapped was successful
  if(max(genesets_table_final$prop_unmapped_items) > max_prop_unmapped){
    stop(glue("Sets were not filtered for max_prop_unmapped: {max_prop_unmapped.}."))
  }

  # test that the final size of sets is not larger than the initial one
  genesets_table_final_sets_size <- genesets_table_final[,c("annotation_id", "size", "size_initial_set")]
  if(!all(genesets_table_final_sets_size$size_initial_set >= genesets_table_final_sets_size$size)){
    stop("Final size of sets cannot be larger than the initial one.")
  }

  # Check expected number of final sets
  n_unique_gene_sets <- length(unique(msig_table$gs_id))
  n_final_gene_sets <- length(unique(genesets_table_final$annotation_id))
  if(n_final_gene_sets != (n_unique_gene_sets-length(sets_removed))){
    stop(glue("The number of final sets: {n_final_gene_sets} does not match expectations of initial sets: {n_unique_gene_sets} - sets removed {sets_removed}"))
  }

  # test name and ids of final sets should be the same as initial - check that any mergin didn't mess things up
  msig_table_unique_sets <- msig_table[,c("gs_id", "gs_name")] %>% unique()
  msig_table_unique_sets <- msig_table_unique_sets[!(msig_table_unique_sets$gs_id %in% sets_removed),]

  genesets_table_final_unique_ids <- genesets_table_final[,c("annotation_id", "annotation_name")] %>% unique()
  stopifnot(nrow(genesets_table_final_unique_ids) == nrow(msig_table_unique_sets))

  order_sets <- match(msig_table_unique_sets$gs_id, genesets_table_final_unique_ids$annotation_id)
  genesets_table_final_unique_ids <- genesets_table_final_unique_ids[order_sets, ]

  id_compare <- all(genesets_table_final_unique_ids$annotation_id != msig_table_unique_sets$gs_id)
  name_compare <- all(genesets_table_final_unique_ids$annotation_name != msig_table_unique_sets$gs_name)
  if(id_compare & name_compare){
    stop(glue("Annotation names or ids don't match with initial MsigDB table. Name: {name_compare}; Ids: {id_compare}"))
  }

  info(EnrichmentLogger(), "Check correspondance of final table and gene sets.")
  # check correspondance of genesets table and items table
  unique_genesets_table <- unique(genesets_table_final$annotation_id)
  unique_gene_sets_genes <- unique(genesets_items_final$annotation_id)

  stopifnot(all(unique_genesets_table %in% unique_gene_sets_genes))
  stopifnot(all(unique_gene_sets_genes %in% unique_genesets_table))
  stopifnot(length(unique_gene_sets_genes) == length(unique_genesets_table))
}


#' define_genesets_to_remove
#'
#' @param genesets_table Table with gene sets and their information
#' @param max_prop_unmapped unmapped proportion to use for filtering gene sets
#'
#' @keywords internal
define_genesets_to_remove <- function(genesets_table, max_prop_unmapped){
  remove_gene_sets <- genesets_table$gs_id[genesets_table$PropUnmapped > max_prop_unmapped]
  return(remove_gene_sets)
}

#' rename_geneset_table
#'
#' @param genesets_table Table with gene sets and their infos
#'
#' @importFrom dplyr rename
#' @export
rename_genesets_table_columns <- function(genesets_table){
  genesets_table <- genesets_table |>
    dplyr::rename(annotation_name = gs_name,
                  annotation_id = gs_id,
                  organism = gs_org,
                  annotation_category = gs_cat,
                  annotation_subcategory = gs_subcat,
                  annotation_exact_source = gs_exact_source,
                  linkout = gs_url,
                  annotation_description = gs_description,
                  annotation_description_full = gs_description_full,
                  size_initial_set = itemsNInitial,
                  size = itemsNUniprot,
                  prop_unmapped_items = PropUnmapped
                  )
  return(genesets_table)
}

#' rename_genesets_items
#'
#' @param genesets_items Table with mapping between gene sets and their items
#'
#' @importFrom dplyr rename
#' @export
rename_genesets_items_columns <- function(genesets_items){
  genesets_items <- genesets_items |> dplyr::rename(annotation_id = gs_id)
  return(genesets_items)
}

#' rename_genes_protein_mapping
#'
#' @param genes_protein_mapping Table with gene sets, gene names, and uniprot id mapping
#'
#' @importFrom dplyr rename
#' @export
rename_genes_protein_mapping_columns <- function(genes_protein_mapping){
  genes_protein_mapping <- genes_protein_mapping |>
    dplyr::rename(annotation_id = gs_id,
                  GeneName = items,
                  UniprotId = UNIPROTID,
                  item_in_initial_set = itemsInitial,
                  item_in_uniprot_set = itemsUniprot,
                  version_uniprot = version
    )
  return(genes_protein_mapping)
}


#' generate_msig_sets
#'
#' @param xml_filepath file path to msigDB xml file
#'
#' @import dplyr
#' @import xml2
#' @export
msigdb_xml_to_table <- function(xml_filepath){
  # Import the MSigDB XML file and extract fields of interest
  info(EnrichmentLogger(), glue("MsigDB xml file exist: {file.exists(xml_filepath)}"))
  mdb_doc <- read_xml(xml_filepath)

  # Code from https://github.com/igordot/msigdbr/blob/main/data-raw/msigdbr-prepare.R
  # MEMBERS_SYMBOLIZED
  # MEMBERS_MAPPING: ensemble IDs
  # MEMBERS_EZID: entrez ids

  mdb_gs_ns <- xml_find_all(mdb_doc, xpath = ".//GENESET")
  mdb_tbl <-
    tibble(
      gs_name = xml_attr(mdb_gs_ns, attr = "STANDARD_NAME"),
      gs_id = xml_attr(mdb_gs_ns, attr = "SYSTEMATIC_NAME"),
      gs_org = xml_attr(mdb_gs_ns, attr = "ORGANISM"),
      gs_cat = xml_attr(mdb_gs_ns, attr = "CATEGORY_CODE"),
      gs_subcat = xml_attr(mdb_gs_ns, attr = "SUB_CATEGORY_CODE"),
      gs_exact_source = xml_attr(mdb_gs_ns, attr = "EXACT_SOURCE"),
      gs_url = xml_attr(mdb_gs_ns, attr = "EXTERNAL_DETAILS_URL"),
      gs_description = xml_attr(mdb_gs_ns, attr = "DESCRIPTION_BRIEF"),
      gs_description_full = xml_attr(mdb_gs_ns, attr = "DESCRIPTION_FULL"),
      gs_members = xml_attr(mdb_gs_ns, attr = "MEMBERS_SYMBOLIZED")
    ) %>%
    dplyr::filter(gs_cat != "ARCHIVED")

  return(mdb_tbl)

}

#' create_gene_set_table
#'
#' @param msigdb_table tibble with data imported with `msigdb_xml_to_table()`
#'
#' @import dplyr
#' @export
create_gene_set_table <- function(msigdb_table){
  # Create a table for gene sets
  msigdbr_genesets <- msigdb_table %>%
    select(!gs_members) %>%
    distinct() %>%
    arrange(gs_name, gs_id)
  return(msigdbr_genesets)
}

#' map_items_to_uniprot
#'
#' @param msigdb_table tibble with data imported with `msigdb_xml_to_table()`
#' @param uniprot_table_path File path to table downloaded from uniprot
#'
#' @import dplyr
#' @import data.table
#' @export
map_items_to_uniprot <- function(msigdb_table, uniprot_table_path){
  uniprot_table <- fread(uniprot_table_path)
  name_version <- basename(uniprot_table_path)

  # Create a table for genes in a tidy/long format (one gene per row)
  geneset_genes <- dplyr::select(msigdb_table, gs_id, gs_members)
  geneset_genes <- dplyr::mutate(geneset_genes, items = strsplit(gs_members, ",", fixed = TRUE))
  geneset_genes <- geneset_genes |> select(gs_id, items) |> tidyr::unnest(cols = items)

  mapping_genes_to_uniprot <- map_symbols_to_uniprotids(items = unique(geneset_genes$items), uniprot_table)
  mapped <- mapping_genes_to_uniprot$mapped_genes

  stopifnot(table(table(mapped$UNIPROTID)) == nrow(mapped))

  geneset_genes <- geneset_genes |> left_join(mapped, by = c("items" = "GeneName"))
  geneset_genes <- geneset_genes |> mutate(itemsInitial = !is.na(items), itemsUniprot = !is.na(UNIPROTID))
  geneset_genes$version <- name_version
  return(geneset_genes)
}

#' compute_unmapped_prop
#'
#' @param geneset_genes table with gene sets and proteins/genes mapped to them
#' @import dplyr
#' @export
compute_unmapped_prop <- function(geneset_genes){
  msigdb_geneset_members <- geneset_genes |>
    dplyr::group_by(gs_id) |>
    dplyr::summarise(itemsNInitial = sum(itemsInitial),
              itemsNUniprot = sum(itemsUniprot))
  msigdb_geneset_members$PropUnmapped <- (msigdb_geneset_members$itemsNInitial - msigdb_geneset_members$itemsNUniprot)/msigdb_geneset_members$itemsNInitial
  return(msigdb_geneset_members)
}


#' map_symbols_to_uniprotids
#'
#' @description Takes a list of genes symbols, the UNIPROT table containing both the "Entry" (uniprot id) and the 'Gene Names (primary)' columns
#' and map gene symbols to uniprotIDs.

#' @param items list of genes
#' @param uniprot_table Table downloaded from uniprot
#'
#' @import dplyr
#' @export
map_symbols_to_uniprotids <- function(items, uniprot_table){

  # case of gene sets with no genes
  if(length(items) == 0){
    return(list(mapped_genes = NULL, unmapped_genes = NULL))
  }

  # otherwise, check the uniprot table and filter for the genes that we need
  if(!(all(c("Entry","Reviewed","Gene Names (primary)") %in% colnames(uniprot_table)))){
    stop("Uniprot table must contain the following columns: 'Entry', 'Reviewed', 'Gene Names (primary)'")
  }

  # uniprot_table <- uniprot_table[uniprot_table$Reviewed %in% "reviewed", ]
  mapping_genes <- uniprot_table[uniprot_table$`Gene Names (primary)` %in% items,
                                 c("Entry","Reviewed","Gene Names (primary)")]
  mapping_genes <- mapping_genes %>% dplyr::rename(GeneName = `Gene Names (primary)`, UNIPROTID = Entry)

  # case when no genes are mapped
  if(nrow(mapping_genes) == 0){
    stop("No genes were mapped to uniprot ids.")
  }

  # keep only reviewed uniprot proteins to remove multimappings
  unique_genes <- unique(mapping_genes$GeneName)
  filter_mappings <- lapply(unique_genes, function(gene) extract_one_entry_uniprot_kb(mapping_genes[mapping_genes$GeneName == gene,]))
  filter_mappings_df <- bind_rows(filter_mappings)

  # some genes might be filtered out as they don't belong to reviewed proteins so we get more unmapped genes
  unmapped <- items[!(items %in% mapping_genes$GeneName)]

  freq_mapped_proteins <- table(filter_mappings$GeneName)
  multimapping_genes <- names(freq_mapped_proteins)[freq_mapped_proteins > 1]
  multimapping_genes_df <- NULL
  if(length(multimapping_genes) > 0){
    time_rep <- freq_mapped_proteins[names(freq_mapped_proteins) %in% multimapping_genes]
    multimapping_genes_df <- data.frame(GeneName = multimapping_genes, Time = time_rep)
  }

  # do we keep both proteins that are mapped by the same gene?
  if(nrow(filter_mappings_df) > length(items)){
    info(EnrichmentLogger(), "The number of mapped genes is larger than the number of initial items.")
  }

  return(list(mapped_genes = filter_mappings_df,
              multimapping_genes=multimapping_genes_df))

}


#' extract_one_entry_uniprot_kb
#' @keywords internal
extract_one_entry_uniprot_kb <- function(map_one_gene){
  map_one <- map_one_gene
  if(nrow(map_one_gene) > 1 & ("reviewed" %in% map_one_gene$Reviewed)){
    map_one <- map_one_gene[map_one_gene$Reviewed %in% "reviewed",]
  }
  return(map_one)
}

