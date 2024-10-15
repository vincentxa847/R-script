library(dplyr)
library(openxlsx)

# -----------------------------------------------------------------------------
#### Prepare gene list ####
# -----------------------------------------------------------------------------
# Function to read gene lists from text files and handle errors
read_gene_list <- function(file_path) {
  if (file.exists(file_path)) {
    if (grepl("HPOList", file_path)) {
      # For files in "HPOList", read with header and select "name" column (this is directly download from HPO website)
      gene_list <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)$name
      gene_list <- trimws(gene_list)
    } else {
      # For other files, read without a header and select the first column
      gene_list <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
    }
    return(gene_list)  
  } else {
    stop(paste("File not found:", file_path))  
  }
}

# Define paths for the gene lists
# TODO: REMAIN TO UPDATE
gene_paths <- list(
  top_candidate_genes = "../GENE_LIST/TOP_candidate_1006.txt",
  candidate_genes = "../GENE_LIST/candidate_1006.txt",
  top_candidate_related_genes = "../GENE_LIST/TOP_candidate_STRING_1006.txt",
  # SCGS
  cilium = "../GENE_LIST/cilium_components.txt", 
  # KEGG
  ECM_interaction = "../GENE_LIST/KEGG/hsa04512_ECM_receptor_interaction.txt",
  hsa04310_Wnt = "../GENE_LIST/KEGG/hsa04310_Wnt.txt",
  hsa04330_Notch = "../GENE_LIST/KEGG/hsa04330_Notch.txt",
  hsa04340_hedgedog = "../GENE_LIST/KEGG/hsa04340_hedgedog.txt",
  hsa04020_calcium_signaling = "../GENE_LIST/KEGG/hsa04020_calcium_signaling.txt",
  hsa04350_TGFbeta = "../GENE_LIST/KEGG/hsa04350_TGFbata.txt",
  N01453_BMP = "../GENE_LIST/KEGG/N01453_BMP.txt",
  # HPO
  HP0006695_ECD = "../GENE_LIST/HPOList/genes_for_HP_0006695",
  HP0001627_CHD = "../GENE_LIST/HPOList/genes_for_HP_0001627.txt",
  # GO
  protein_folding = "../GENE_LIST/GO/GO0140662_ATPdependent protein folding chaperone.txt"
)

# Read in gene list
gene_lists <- lapply(gene_paths, read_gene_list)

# -----------------------------------------------------------------------------
#### Output a list of gene for IPA #### 
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# TODO: Investigate how to use IPA (Ingenuity Pathway Analysis) for variant analysis
#       and determine what additional information should be included.
# -----------------------------------------------------------------------------

# Function to extract gene symbols from variants and save them to a text file
# If more than one variant exists in a gene, multiple gene symbols may be written.
# These duplicates will be handled (removed) during IPA processing.
writeout_list <- function(table, gene_list, output_name) {
  output <- table %>% filter(Gene_refgene %in% gene_list)
  output <- output$Gene_refgene 
  
  file_path <- paste0("../IPA/", output_name, ".txt")
  
  write.table(output, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Example usage : 
out <- writeout_list(D25029,gene_lists$HP0006695_ECD, "D25029_MAF0.01_HP0006695")  

# -----------------------------------------------------------------------------
#### Summary list and EXCEL output of variants in different genelist ####
# -----------------------------------------------------------------------------
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

# Function to summarize the data by Gene and Variant Type, which is a integration of 05.Summarize.R
summarize_by_gene_variant <- function(df) {
  # Group by 'Gene_refgene' and 'Func_refgene', then summarize counts
  summary_table <- df %>%
    group_by(Gene_refgene, Func_refgene) %>%
    summarize(count = n(), .groups = "drop")
  
  # Pivot table to a wider format, fill missing values with 0
  wide_table <- summary_table %>%
    tidyr::pivot_wider(names_from = Func_refgene, values_from = count, values_fill = 0)
  
  # Ensure the required columns ('exonic', 'UTR3', 'UTR5', 'intronic', 'upstream', 'downstream', 'intergenic') are present
  required_columns <- c("exonic", "UTR3", "UTR5", "intronic", "upstream", "downstream", "intergenic")
  # Add missing columns with default value 0
  for (col in required_columns) {
    if (!col %in% colnames(wide_table)) {
      wide_table[[col]] <- 0
    }
  }
  # Reorder columns to ensure the required columns are at the front
  wide_table <- wide_table %>%
    dplyr::select(Gene_refgene, all_of(required_columns), everything())
  
  return(wide_table)
}

# Function to include variants in gene list and produce EXCEL output
gene_list <- function(input, gene_list, output_name) {
  
  if (!"Gene_refgene" %in% colnames(input)) {
    stop("The input data does not contain a column named 'Gene_refgene'.")
  }

  # Clean the data before processing
  input <- clean_illegal_strings(input)
  
  # Define a list to hold all gene_list_variant
  gene_list_variant <- list()
  
  # 1.top_candidate_genes
  gene_list_variant$top_candidate_genes <- input %>% filter(Gene_refgene %in% gene_list$top_candidate_genes)
  # 2.top_candidate_related_genes
  gene_list_variant$top_candidate_related_genes <- input %>% filter(Gene_refgene %in% gene_list$top_candidate_related_genes)
  # 3.candidate_genes
  gene_list_variant$candidate_genes <- input %>% filter(Gene_refgene %in% gene_list$candidate_genes)
  # 4.cilium
  gene_list_variant$cilium <- input %>% filter(Gene_refgene %in% gene_list$cilium)
  # 5.ECM_interaction
  gene_list_variant$ECM_interaction <- input %>% filter(Gene_refgene %in% gene_list$ECM_interaction)
  # 6.hsa04310_Wnt
  gene_list_variant$hsa04310_Wnt <- input %>% filter(Gene_refgene %in% gene_list$hsa04310_Wnt)
  # 7.hsa04330_Notch
  gene_list_variant$hsa04330_Notch <- input %>% filter(Gene_refgene %in% gene_list$hsa04330_Notch)
  # 8.hsa04340_hedgedog
  gene_list_variant$hsa04340_hedgedog <- input %>% filter(Gene_refgene %in% gene_list$hsa04340_hedgedog)
  # 9.hsa04020_calcium_signaling
  gene_list_variant$hsa04020_calcium_signaling <- input %>% filter(Gene_refgene %in% gene_list$hsa04020_calcium_signaling)
  # 10.hsa04350_TGFbata
  gene_list_variant$hsa04350_TGFbeta <- input %>% filter(Gene_refgene %in% gene_list$hsa04350_TGFbeta)
  # 11.N01453_BMP
  gene_list_variant$N01453_BMP <- input %>% filter(Gene_refgene %in% gene_list$N01453_BMP)
  # 12.HP0006695_ECD
  gene_list_variant$HP0006695_ECD <- input %>% filter(Gene_refgene %in% gene_list$HP0006695_ECD)
  # 13.HP0001627_CHD
  gene_list_variant$HP0001627_CHD <- input %>% filter(Gene_refgene %in% gene_list$HP0001627_CHD)
  # 14. chr21
  gene_list_variant$chr21 <- input %>% filter(Chr == "21")
  
  ## Prepare genesymbol worksheet and Sort each group by CADD_phred
  # Create a list to hold Gene_refgene columns for all groups
  gene_refgene_list <- list()
  # Loop through each variant group and extract Gene_refgene
  for (group_name in names(gene_list_variant)) {
    
    # Sort each group by CADD_phred
    if ("CADD_phred" %in% colnames(gene_list_variant[[group_name]])) {
      gene_list_variant[[group_name]] <- gene_list_variant[[group_name]] %>%
        arrange(desc(CADD_phred))  # Sort in descending order by CADD_phred
    }
    
    gene_column <- gene_list_variant[[group_name]]$Gene_refgene
    
    # Fill with empty string for unequal lengths
    gene_refgene_list[[group_name]] <- c(gene_column, rep("", max(0, max(sapply(gene_list_variant, nrow)) - length(gene_column))))
  }
  # Convert the list to a data frame
  gene_refgene_df <- as.data.frame(gene_refgene_list)

  ## Prepare summary worksheet
  summarized_results <- lapply(gene_list_variant, summarize_by_gene_variant)
  names(summarized_results) <- names(gene_list_variant)
  
  ## EXCEL output
  wb <- createWorkbook() # Create a workbook for the Excel output
  # write genesymbol worksheet
  addWorksheet(wb, "GENESYMBOL OF GENELIST")
  writeData(wb, "GENESYMBOL OF GENELIST", gene_refgene_df)
  # write summary worksheet
  start_row <- 1
  addWorksheet(wb, "Summary")
  for (table in names(summarized_results)) {
    writeData(wb, "Summary", paste(table), startRow = start_row, startCol = 1) # Write table name
    start_row <- start_row + 1
    writeData(wb, "Summary", summarized_results[[table]], startRow = start_row, startCol = 1) # Write data frame
    start_row <- start_row + nrow(summarized_results[[table]]) + 2 # Move start_row down for the next table
  }
  
  # Create individual worksheets for each variant group
  for (group_name in names(gene_list_variant)) {
    addWorksheet(wb, group_name)  # Create a new worksheet for each group
    writeData(wb, group_name, gene_list_variant[[group_name]])  # Write the group data to the worksheet
  }
  
  # Save the workbook to the specified file
  saveWorkbook(wb, output_name, overwrite = TRUE)
  
  # Return the list of grouped variants
  return(gene_list_variant)
}

# -----------------------------------------------------------------------------
#### Usage ####
# -----------------------------------------------------------------------------
D25007_GeneList <- gene_list(D25007_MAF0.01,gene_lists,"../nonDS-ECD/D25007/D25007_MAF0.01.gene_lists.xlsx")
D25029_GeneList <- gene_list(D25029_MAF0.01,gene_lists,"../DS-ECD/D25029/D25029_MAF0.01.gene_lists.xlsx")
D25046_GeneList <- gene_list(D25046_MAF0.01,gene_lists,"../DS-ECD/D25046/D25046_MAF0.01.gene_lists.xlsx")

D25165_GeneList <- gene_list(D25165_MAF0.01,gene_lists,"../DS-nonECD/D25165/D25165_MAF0.01.gene_lists.xlsx")
D25168_GeneList <- gene_list(D25168_MAF0.01,gene_lists,"../DS-nonECD/D25168/D25168_MAF0.01.gene_lists.xlsx")

