library(dplyr) 
library(tidyr)
library(openxlsx)
library(ggplot2)
library(ggVennDiagram)

# -------------------------------------------------------------------
# Compare and visualize the shared records
# -------------------------------------------------------------------
# Create the ID column for each sample
# DS-ECD
D25029_MAF0.01 <- D25029_MAF0.01 %>%
  mutate(ID = paste(Gene.refgene, Chr, Start, End, Ref, Alt, sep = "_"))
D25046_MAF0.01 <- D25046_MAF0.01 %>%
  mutate(ID = paste(Gene.refgene,Chr, Start, End, Ref, Alt, sep = "_"))
# DS-nonECD
D25163_MAF0.01 <- D25163_MAF0.01 %>%
  mutate(ID = paste(Gene.refgene,Chr, Start, End, Ref, Alt, sep = "_"))
D25165_MAF0.01 <- D25165_MAF0.01 %>%
  mutate(ID = paste(Gene.refgene,Chr, Start, End, Ref, Alt, sep = "_"))
D25168_MAF0.01 <- D25168_MAF0.01 %>%
  mutate(ID = paste(Gene.refgene,Chr, Start, End, Ref, Alt, sep = "_"))
# nonDS-ECD
D25007_MAF0.01 <- D25007_MAF0.01 %>%
  mutate(ID = paste(Gene.refgene,Chr, Start, End, Ref, Alt, sep = "_"))

# Variant Id within group to plot
venn_data <- list(
  DS_ECD = c(D25029_MAF0.01$ID, D25046_MAF0.01$ID),
  DS_nonECD = c(D25163_MAF0.01$ID, D25165_MAF0.01$ID, D25168_MAF0.01$ID),
  nonDS_ECD = c(D25007_MAF0.01$ID)
)

# Venn diagram
ggVennDiagram(venn_data, color = "black", 
              lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_x_continuous(expand = expansion(mult = .2))


### Extract shared within groups ####
### Shared variants in three groups (11993)
shared_variants_threeID <- Reduce(intersect, venn_data)
shared_variants_threeID <- trimws(shared_variants_threeID)
shared_variants_three <- D25007_MAF0.01 %>% filter(ID %in% shared_variants_three)

### Shared variants in DS_ECD and nonDS_ECD and exclude 3 overlap (4743)
shared_variants_DSECD_nonDSECDID <- intersect(venn_data$DS_ECD, venn_data$nonDS_ECD)
shared_variants_DSECD_nonDSECDID <- setdiff(shared_variants_DSECD_nonDSECD,shared_variants_threeID)
shared_variants_DSECD_nonDSECD <- rbind(
  D25029_MAF0.01 %>% filter(ID %in% shared_variants_DSECD_nonDSECDID),
  D25046_MAF0.01 %>% filter(ID %in% shared_variants_DSECD_nonDSECDID)
) %>% distinct(ID, .keep_all = TRUE)

### Shared variants in DS_ECD and DS_nonECD and exclude 3 overlap (26471)
shared_variants_DSECD_DSnonECDID <- intersect(venn_data$DS_ECD, venn_data$DS_nonECD)
shared_variants_DSECD_DSnonECDID <- setdiff(shared_variants_DSECD_DSnonECDID,shared_variants_threeID)
shared_variants_DSECD_DSnonECD <- rbind(
  D25029_MAF0.01 %>% filter(ID %in% shared_variants_DSECD_DSnonECDID),
  D25046_MAF0.01 %>% filter(ID %in% shared_variants_DSECD_DSnonECDID)
) %>% distinct(ID, .keep_all = TRUE)

# -------------------------------------------------------------------
# Function to export shared variants
# -------------------------------------------------------------------
groupandExport  <- function(tables, output_name) {

  wb <- createWorkbook() 
  summary_tables_list <- list()
  for (name_ in names(tables)) {
    
    # Create worksheet for each table
    addWorksheet(wb, name_)
    
    ## summary of each table, add at the top of each table
    summary_table <- tables[[name_]] %>%
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
    
    # Save the summary table for the "SUMMARY" worksheet
    summary_tables_list[[name_]] <- wide_table
    
    ## Write the full variant table
    writeData(wb, name_, tables[[name_]])
  }
  addWorksheet(wb, "SUMMARY")
  
  # Combine all summary tables and order by exonic, UTR3, UTR5 in descending order
  combined_summary_table <- bind_rows(summary_tables_list, .id = "Source") %>%
    arrange(desc(exonic), desc(UTR3), desc(UTR5))
  
  writeData(wb, "SUMMARY", combined_summary_table, startRow = 1, startCol = 1)
  
  saveWorkbook(wb, output_name, overwrite = TRUE)
}

groupandExport(list(
  shared_variants_three = shared_variants_three,
  shared_variants_DSECD_DSnonECD = shared_variants_DSECD_DSnonECD,
  shared_variants_DSECD_nonDSECD = shared_variants_DSECD_nonDSECD
), "sharedVariant.xlsx")
