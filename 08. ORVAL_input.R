#### Data for ORVAL Digenic and Oligogenic Analysis ####
library(openxlsx)
library(dplyr)

# Function to clean illegal characters from the data
clean_illegal_strings <- function(df) {
  df_clean <- df
  
  # Loop through each column
  for (col in colnames(df_clean)) {
    if (is.character(df_clean[[col]])) {
      # Remove non-printable characters and invalid UTF-8 strings
      df_clean[[col]] <- iconv(df_clean[[col]], from = "UTF-8", to = "UTF-8", sub = "")
      df_clean[[col]] <- gsub("[[:cntrl:]]", "", df_clean[[col]])  # Remove control characters
    }
  }
  
  return(df_clean)
}

# Function to extract and transform data for ORVAL analysis
ORVAL_variantList <- function(orval_datasets, output_file) {
  # Create a new Excel workbook
  wb <- createWorkbook()
  
  # Iterate over each dataset in orval_datasets
  for (dataset_name in names(orval_datasets)) {
    data <- orval_datasets[[dataset_name]]
    
    # Extract relevant columns, filter, and modify Genotype values
    # ORVAL requires five columns: Chr, Start, Ref, Alt, and Zygosity (Genotype).
    data_modified <- data %>%
      select(Gene_refgene, Func_refgene, Chr, Start, Ref, Alt, Genotype) %>%          # Select required columns
      filter(Ref != 0 & Alt != 0) %>%                     # Filter out rows with invalid values
      mutate(Genotype = case_when(
        Genotype == "het" ~ "Heterozygous",               # Convert 'het' to 'Heterozygous'
        Genotype == "hom" ~ "Homozygous",                 # Convert 'hom' to 'Homozygous'
        Genotype == "hem" ~ "Homozygous",                 # Convert 'hem' to 'Homozygous'
        TRUE ~ Genotype                                   # Leave other values unchanged
      ))
    
    # Clean the data before writing
    data_modified <- clean_illegal_strings(data_modified)
    
    # Create a new worksheet with the dataset name
    addWorksheet(wb, dataset_name)  # Use the dataset name as the worksheet name
    writeData(wb, dataset_name, data_modified)
  }
  
  # Save the workbook to the specified file
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  # Return a success message for clarity
  cat("Data written to", output_file, "\n")
}

#### Apply the Function to Multiple dataset and incorporate into an EXCEL workbook ####
# List of dataset with corresponding names
orval_datasets <- list(
  "D25029.ORVAL.exonic.MAF0.01" = D25029_variantgroup$exonic,
  "D25029.ORVAL.intronic.MAF0.01" = D25029_variantgroup$intronic,
  "D25029.ORVAL.UTR.MAF0.01" = D25029_variantgroup$UTR,
  "D25029.ORVAL.promoter.MAF0.01" = D25029_variantgroup$promoter,
  "D25029.ORVAL.chr21.MAF0.01" = D25029_GeneList$chr21
)

# Write specified dataset to an Excel workbook
ORVAL_writeput <- ORVAL_variantList(orval_datasets, "../ORVAL/D25029.ORVAL.xlsx")



# -----------------------------------------------------------------------------
#### Analysis of chr21 variants with exonic variants in other chromosome ####
# -----------------------------------------------------------------------------
# D25029: ≥ Disease-causing_99.9%-zone, total 3918 entries
# D25046: ≥ Disease-causing_99.9%-zone, total 2655 entries
ORVAL_D25046Chr21Exonic <- read.csv(file = "../ORVAL_RESULT/D25046Chr21Exonic.tsv",sep = "\t") 

# Filter Variant Pair with variant on chr21
# D25029: 46 entries, D25046: 10 entries
ORVAL_D25046Chr21Exonic_ONLYW21 <- ORVAL_D25046Chr21Exonic[
  grepl("^21:", ORVAL_D25046Chr21Exonic$GeneA_Alleles) |
    grepl("^21:", ORVAL_D25046Chr21Exonic$GeneB_Alleles), 
]


# CHD filter (ECD filter result in no data)
ORVAL_D25046Chr21Exonic_ONLYW21_CHD <- ORVAL_D25046Chr21Exonic_ONLYW21[
    (ORVAL_D25046Chr21Exonic_ONLYW21$Gene_A %in% gene_lists$HP0001627_CHD |
       ORVAL_D25046Chr21Exonic_ONLYW21$Gene_B %in% gene_lists$HP0001627_CHD),
]

# Export the filtered data to a TSV file
write.table(ORVAL_D25046Chr21Exonic_ONLYW21_CHD, 
            file = "ORVAL_D25046Chr21Exonic_ONLYW21_CHD.tsv", 
            sep = "\t",        # Use tab as the separator
            row.names = FALSE,  # Do not write row names
            quote = FALSE)      # Do not include quotes around the data
