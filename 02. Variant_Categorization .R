library(dplyr)
library(openxlsx)

group_variants_and_save <- function(data,output_name) {
  
  # Define a list to hold all groups of variants
  variant_groups <- list()
  
  ## Grouping data into 14 distinct groups
  # 1. CADD score > 30
  variant_groups$CADD30 <- subset(data, round(as.numeric(CADD_phred), 2) >= 30)
  
  # 2. CADD score > 15
  variant_groups$CADD15 <- subset(data, round(as.numeric(CADD_phred), 2) >= 15)
  
  # 3. Missense variants with CADD > 15
  variant_groups$CADD15_missense <- subset(data, round(as.numeric(CADD_phred), 2) >= 15 & ExonicFunc_refgene == "nonsynonymous SNV")
  
  # 4. All missense variants
  variant_groups$missense <- subset(data, ExonicFunc_refgene == "nonsynonymous SNV")
  
  # 5. Synonymous variants with CADD > 15
  variant_groups$synonymous_CADD15 <- subset(data, ExonicFunc_refgene == "synonymous SNV" & round(as.numeric(CADD_phred), 2) >= 15)
  
  # 6. All synonymous variants
  variant_groups$synonymous <- subset(data, ExonicFunc_refgene == "synonymous SNV")
  
  # 7. All exonic variants
  variant_groups$exonic <- subset(data, Func_refgene == "exonic")
  
  # 8. All variants in UTR regions
  variant_groups$UTR <- subset(data, Func_refgene == "UTR3" | Func_refgene == "UTR5")
  
  # 9. All intronic variants
  variant_groups$intronic <- subset(data, Func_refgene == "intronic")
  
  # 10. Variants in promoter regions (2kb upstream and downstream of transcriptional start site)
  variant_groups$promoter <- subset(data, Func_refgene == "upstream" | Func_refgene == "downstream")
  
  ## Prepare for conserved regions
  # Replace "." with NA in both phyloP and phastCons columns and convert to numeric
  # some data with NA not ".", so add additional step to replaces any character that is not a digit or a period with NA 
  data$phyloP470way_mammalian_rankscore <- as.numeric(na_if(gsub("[^0-9.]", NA, data$phyloP470way_mammalian_rankscore), "."))
  data$phastCons470way_mammalian_rankscore <- as.numeric(na_if(gsub("[^0-9.]", NA, data$phastCons470way_mammalian_rankscore), "."))
  # Calculate thresholds (75th percentile) for phyloP and phastCons scores
  phyloP_threshold <- quantile(data$phyloP470way_mammalian_rankscore, 0.75, na.rm = TRUE)
  phastCons_threshold <- quantile(data$phastCons470way_mammalian_rankscore, 0.75, na.rm = TRUE)
  
  # 11. Variants in conserved regions
  variant_groups$conserved <- subset(data, phyloP470way_mammalian >= phyloP_threshold & phastCons470way_mammalian_rankscore >= phastCons_threshold)
  
  # 12. Variants in non-coding RNAs (ncRNAs)
  variant_groups$ncRNA <- data[grep("ncRNA", data$Func_refgene),]
  
  # 13. variants in double-elite enhancers (skip now)
  # pass
  
  # 14. loss of function (unsure the criteria, skip now)
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
    
    gene_column <- variant_groups[[group_name]]$Gene_refgene
    
    # Fill with empty string for unequal lengths
    gene_refgene_list[[group_name]] <- c(gene_column, rep("", max(0, max(sapply(variant_groups, nrow)) - length(gene_column))))
  }
  # Convert the list to a data frame
  gene_refgene_df <- as.data.frame(gene_refgene_list)
  
  
  ## EXCEL output
  wb <- createWorkbook() # Create a workbook for the Excel output
  addWorksheet(wb, "GENESYMBOL OF DIFFERENT GROUPS")
  writeData(wb, "GENESYMBOL OF DIFFERENT GROUPS", gene_refgene_df)
  
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


D25029_variantgroup <- group_variants_and_save(D25029,"D25029_variantgroup.xlsx")
D25046_variantgroup <- group_variants_and_save(D25046,"D25046_variantgroup.xlsx")

