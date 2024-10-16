library(dplyr)
library(openxlsx)

# -----------------------------------------------------------------------------
#### Summary list and EXCEL output of variants of different catalog ####
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
  # Group by 'Gene.refgene' and 'Func.refgene', then summarize counts
  summary_table <- df %>%
    group_by(Gene.refgene, Func.refgene) %>%
    summarize(count = n(), .groups = "drop")
  
  # Pivot table to a wider format, fill missing values with 0
  wide_table <- summary_table %>%
    tidyr::pivot_wider(names_from = Func.refgene, values_from = count, values_fill = 0)
  
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
    dplyr::select(Gene.refgene, all_of(required_columns), everything())
  
  return(wide_table)
}

# Function to group variants and produce EXCEL output
group_variants_and_save <- function(data,output_name) {
  
  # Clean the data before processing
  data <- clean_illegal_strings(data)
  
  # Handle non-numeric characters for calculation
  data$CADD_phred <- as.numeric(gsub("[^0-9.]", NA, data$CADD_phred))
  data$phyloP470way_mammalian_rankscore <- as.numeric(gsub("[^0-9.]", NA, data$phyloP470way_mammalian_rankscore))
  data$phastCons470way_mammalian_rankscore <- as.numeric(gsub("[^0-9.]", NA, data$phastCons470way_mammalian_rankscore))
  
  # Define a list to hold all groups of variants
  variant_groups <- list()
  
  ## Grouping data into 14 distinct groups
  # 1. Loss of Function (STOP GAIN or SPLICING or FRAMESGIFT)
  variant_groups$LOF <- subset(data, Consequence == "stop_gained" | Consequence =="splice_acceptor_variant" | Consequence =="splice_donor_variant" | Consequence == "frameshift_variant")
  
  # 2. CADD score > 30
  variant_groups$CADD30 <- subset(data, !is.na(CADD_phred) & round(as.numeric(CADD_phred), 2) >= 30)
  
  # 3. CADD score > 15
  variant_groups$CADD15 <- subset(data, !is.na(CADD_phred) & round(as.numeric(CADD_phred), 2) >= 15)
  
  # 4. Missense variants with CADD > 15
  variant_groups$CADD15_missense <- subset(data, !is.na(CADD_phred) & round(as.numeric(CADD_phred), 2) >= 15 & ExonicFunc.refgene == "nonsynonymous SNV")
  
  # 5. All missense variants
  variant_groups$missense <- subset(data, ExonicFunc.refgene == "nonsynonymous SNV")
  
  # 6. Synonymous variants with CADD > 15
  variant_groups$synonymous_CADD15 <- subset(data, ExonicFunc.refgene == "synonymous SNV" & !is.na(CADD_phred) & round(as.numeric(CADD_phred), 2) >= 15)
  
  # 7. All synonymous variants
  variant_groups$synonymous <- subset(data, ExonicFunc.refgene == "synonymous SNV")
  
  # 8. All exonic variants
  variant_groups$exonic <- subset(data, Func.refgene == "exonic")
  
  # 9. All variants in UTR regions
  variant_groups$UTR <- subset(data, Func.refgene == "UTR3" | Func.refgene == "UTR5")
  
  # 10. All intronic variants
  variant_groups$intronic <- subset(data, Func.refgene == "intronic")
  
  # 11. Variants in promoter regions (2kb upstream and downstream of transcriptional start site)
  variant_groups$promoter <- subset(data, Func.refgene == "upstream" | Func.refgene == "downstream")
  
  # 12. Variants in conserved regions
  # Calculate thresholds (75th percentile) for phyloP and phastCons scores
  phyloP_threshold <- quantile(data$phyloP470way_mammalian_rankscore, 0.75, na.rm = TRUE)
  phastCons_threshold <- quantile(data$phastCons470way_mammalian_rankscore, 0.75, na.rm = TRUE)
  variant_groups$conserved <- subset(data, phyloP470way_mammalian >= phyloP_threshold & phastCons470way_mammalian_rankscore >= phastCons_threshold)
  
  # 13. Variants in non-coding RNAs (ncRNAs)
  variant_groups$ncRNA <- data[grep("ncRNA", data$Func.refgene),]
  
  # 14. variants in double-elite enhancers (skip now)
  # pass
  
  ## Prepare genesymbol worksheet and Sort each group by CADD_phred
  # Create a list to hold Gene_refgene columns for all groups
  gene_refgene_list <- list()
  # Loop through each variant group and extract Gene_refgene
  for (group_name in names(variant_groups)) {
    
    # Sort each group by CADD_phred
    if ("CADD_phred" %in% colnames(variant_groups[[group_name]])) {
      variant_groups[[group_name]] <- variant_groups[[group_name]] %>%
        arrange(desc(CADD_phred))  # Sort in descending order by CADD_phred
    }
    
    gene_column <- variant_groups[[group_name]]$Gene.refgene
    
    # Fill with empty string for unequal lengths
    gene_refgene_list[[group_name]] <- c(gene_column, rep("", max(0, max(sapply(variant_groups, nrow)) - length(gene_column))))
  }
  # Convert the list to a data frame
  gene_refgene_df <- as.data.frame(gene_refgene_list)
  
  ## Prepare summary worksheet
  summarized_results <- lapply(variant_groups, summarize_by_gene_variant)
  names(summarized_results) <- names(variant_groups)
  
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
  for (group_name in names(variant_groups)) {
    addWorksheet(wb, group_name)  # Create a new worksheet for each group
    writeData(wb, group_name, variant_groups[[group_name]])  # Write the group data to the worksheet
  }
  
  # Save the workbook to the specified file
  saveWorkbook(wb, output_name, overwrite = TRUE)
  
  # Return the list of grouped variants
  return(variant_groups)
}

# -----------------------------------------------------------------------------
#### Usage ####
# -----------------------------------------------------------------------------
D25007_variantgroup <- group_variants_and_save(D25007_MAF0.01,"../nonDS-ECD/D25007/D25007_variantgroup.xlsx")
D25029_variantgroup <- group_variants_and_save(D25029_MAF0.01,"../DS-ECD/D25029/D25029_variantgroup.xlsx")
D25046_variantgroup <- group_variants_and_save(D25046_MAF0.01,"../DS-ECD/D25046/D25046_variantgroup.xlsx")

