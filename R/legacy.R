# legacy code




process_protein_intensities <- function(protein_count_and_intensities, by = "Protein"){
  if (by == "Protein"){
    protein_count_and_intensities = process_protein_intensities_by_protein(protein_count_and_intensities)
  }
  else if(by == "Gene"){
    protein_count_and_intensities = process_protein_intensities_by_gene(protein_count_and_intensities)
  }
}


process_protein_intensities_by_protein <- function(protein_count_and_intensities){
  
  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(conditions)) #first level unwrap
  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(intensityValues)) #second level unwrap
  protein_count_and_intensities = protein_count_and_intensities[,c("ProteinId", "name", "log2NInt_ProteinGroupId")]
  protein_count_and_intensities = protein_count_and_intensities %>% group_by(ProteinId,name)  %>% mutate(replicate = row_number())
  protein_count_and_intensities = protein_count_and_intensities %>%
    group_by(ProteinId) %>%
    unite("ID_type", name, replicate,remove = TRUE) %>%
    spread(ID_type,  log2NInt_ProteinGroupId)
  
  protein_count_and_intensities
}
#

process_protein_intensities_by_gene <- function(protein_count_and_intensities){
  
  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(conditions)) #first level unwrap
  protein_count_and_intensities =  unnest(protein_count_and_intensities,cols = c(intensityValues)) #second level unwrap
  print(colnames(protein_count_and_intensities))
  protein_count_and_intensities = protein_count_and_intensities[,c("ProteinId", "GeneName", "name", "log2NInt_ProteinGroupId")]
  protein_count_and_intensities = protein_count_and_intensities %>% group_by(GeneName, ProteinId,name)  %>% mutate(replicate = row_number())
  protein_count_and_intensities = protein_count_and_intensities %>%
    group_by(GeneName, ProteinId) %>%
    unite("ID_type", name, replicate,remove = TRUE) %>%
    spread(ID_type,  log2NInt_ProteinGroupId)
  
  protein_count_and_intensities
}
#