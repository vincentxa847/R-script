library(dplyr) 

# -------------------------------------------------------------------
# Function to identify overlapping variants between two samples
# Also groups the variants by functional regions
# -------------------------------------------------------------------

overlap_variant <- function(table1, table2) {
  # Columns to match: "Chr", "Start", "End", "Ref", "Alt"
  # Perform a merge to find the common rows
  shared_rows <- merge(table1, table2, by = c("Chr", "Start", "End", "Ref", "Alt"))
  
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
  
  return(variant_groups)
}

# Example usage
tmp_overlap_variant <- overlap_variant(D25165_MAF0.01, D25168_MAF0.01)

# Extract the exonic group and select relevant columns
tmp_overlap_variant <- tmp_overlap_variant$exonic[, 1:13]

# View result
tmp_overlap_variant
