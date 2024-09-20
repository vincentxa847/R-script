#### Data output to ORVAL ####
# VCF of WGS is too large, so use variant list
# ORVAL take FIVE COLUMN (Chr, Start, Ref, Alt and Zygosity)
# By doing so, not need to use in-build filter function in ORVAL 
ORVAL_variantList <- function(data, output_file) {
  # Extract relevant columns and modify Genotype values
  data_modified <- data %>%
    select(Chr, Start, Ref, Alt, Genotype) %>%
    filter(Ref != 0 & Alt != 0) %>%  # Exclude rows where Ref or Alt is 0
    mutate(Genotype = case_when(
      Genotype == "het" ~ "Heterozygous",
      Genotype == "hmo" ~ "Homozygous",
      TRUE ~ Genotype  # Keep other values unchanged, if any
    ))
  
  # Write the modified data to a TSV file
  write.table(data_modified, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

ORVAL_variantList_top_candidate_related_genes_D25007_MAF0.01 = ORVAL_variantList(top_candidate_related_genes_D25007_MAF0.01,"top_candidate_related_genes_D25007_MAF0.01.tsv")

ORVAL_variantList_top_candidate_related_genes_N1675_MAF0.01 = ORVAL_variantList(top_candidate_related_genes_N1675_MAF0.01,"top_candidate_related_genes_N1675_MAF0.01.tsv")
