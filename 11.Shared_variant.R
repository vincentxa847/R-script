library(dplyr) 
library(tidyr)
library(openxlsx)
library(ggplot2)
library(ggVennDiagram)

# -------------------------------------------------------------------
# Function to identify overlapping variants between two samples
# Also groups the variants by functional regions
# -------------------------------------------------------------------

overlap_variant <- function(tables, output_name) {
  # Ensure that tables is a list of data frames
  if (!is.list(tables) || length(tables) < 2) {
    stop("Input should be a list with at least two tables.")
  }
  
  # Columns to match: "Chr", "Start", "End", "Ref", "Alt"
  # Perform a merge to find the common rows across all tables
  shared_rows <- Reduce(function(x, y) merge(x, y, by = c("Gene.refgene","Chr", "Start", "End", "Ref", "Alt")), tables)
  
  # Define a list to hold all groups of variants
  variant_groups <- list()
  
  # Use logical indexing for better performance
  # Group shared variants in "exonic", "UTR regions", "intronic", "promoter regions"
  
  # exonic variants
  variant_groups$exonic <- shared_rows[shared_rows$Func_refgene.x == "exonic", ]
  
  # UTR regions (combining 3'UTR and 5'UTR)
  variant_groups$UTR <- shared_rows[shared_rows$Func_refgene.x %in% c("UTR3", "UTR5"), ]
  
  # intronic variants
  variant_groups$intronic <- shared_rows[shared_rows$Func_refgene.x == "intronic", ]
  
  # promoter regions (upstream or downstream of transcription start site)
  variant_groups$promoter <- shared_rows[shared_rows$Func_refgene.x %in% c("upstream", "downstream"), ]
  
  # EXCEL output
  wb <- createWorkbook()
  for (group_name in names(variant_groups)) {
    addWorksheet(wb, group_name)  
    writeData(wb, group_name, variant_groups[[group_name]])  
  }
  # Save the workbook to the specified file
  saveWorkbook(wb, output_name, overwrite = TRUE)
  
  return(variant_groups)
}

# -------------------------------------------------------------------
# Example usage
# REMAIN TO UPDATE: EXCEL OUTPUT and GROUP-LEVEL COMPARSION
# -------------------------------------------------------------------
DSECD_0.01 <- list(D25029_MAF0.01, D25046_MAF0.01)
DSnonECD_0.01 <- list(D25165_MAF0.01, D25168_MAF0.01)
tmp_overlap_variant <- overlap_variant(DSnonECD_0.01,"../DS-nonECD/DSnonECD_0.01_shareVariant.xlsx")

# -------------------------------------------------------------------
# Compare two tables in R and visualize the shared and unique rows
# need to consider a variant in more then one gene
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


### Extract shared and unique variants
# Find unique rows in D25046_MAF0.01
unique_table2 <- anti_join(D25046_MAF0.01, D25029_MAF0.01, by = "ID")


# Find shared variants of three groups (11993)
shared_variants <- Reduce(intersect, venn_data)
shared_variants <- trimws(shared_variants)
# Using one group can fetch the shared variants of three groups 
shared_variants_data <- D25007_MAF0.01 %>% filter(ID %in% shared_variants)

