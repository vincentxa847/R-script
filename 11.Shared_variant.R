library(dplyr) 
library(openxlsx)

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
  shared_rows <- Reduce(function(x, y) merge(x, y, by = c("Chr", "Start", "End", "Ref", "Alt")), tables)
  
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
tmp_overlap_variant <- overlap_variant(DSECD_0.01,"../DS-ECD/DSECD_0.01_shareVariant.xlsx")
