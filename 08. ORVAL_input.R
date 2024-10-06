#### Data Output to ORVAL ####
# ORVAL requires five columns: Chr, Start, Ref, Alt, and Zygosity.

# This function modifies the data to fit the ORVAL format.
ORVAL_variantList <- function(data, output_file) {
  
  # Extract relevant columns, filter, and modify Genotype values
  data_modified <- data %>%
    select(Chr, Start, Ref, Alt, Genotype) %>%          # Select required columns
    filter(Ref != 0 & Alt != 0) %>%                     # Filter out rows with invalid values
    mutate(Genotype = case_when(
      Genotype == "het" ~ "Heterozygous",               # Convert 'het' to 'Heterozygous'
      Genotype == "hmo" ~ "Homozygous",                 # Convert 'hmo' to 'Homozygous'
      TRUE ~ Genotype                                   # Leave other values unchanged
    ))
  
  # Write the modified data to a TSV file
  write.table(data_modified, file = output_file, sep = "\t", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Return a success message for clarity
  cat("Data written to", output_file, "\n")
}

#### Apply the Function to Multiple Datasets ####
# List of datasets with corresponding output file names
orval_datasets <- list(
  "ORVAL_top_candidate_related_genes_D25029_MAF0.01" = top_candidate_related_genes_D25029_MAF0.01
)

# Apply the ORVAL_variantList function to each dataset and save the result
lapply(names(orval_datasets), function(dataset_name) {
  ORVAL_variantList(orval_datasets[[dataset_name]], paste0(dataset_name, ".tsv"))
})
