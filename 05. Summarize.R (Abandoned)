library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Function to summarize the data by Gene and Variant Type
# -----------------------------------------------------------------------------
summarize_by_gene_variant <- function(df) {
  # Group by 'Gene_refgene' and 'Func_refgene', then summarize counts
  summary_table <- df %>%
    group_by(Gene_refgene, Func_refgene) %>%
    summarize(count = n(), .groups = "drop")
  
  # Pivot table to a wider format, fill missing values with 0
  wide_table <- summary_table %>%
    tidyr::pivot_wider(names_from = Func_refgene, values_from = count, values_fill = 0)
  
  return(wide_table)
}

# -----------------------------------------------------------------------------
# List of datasets to summarize
# TODO: CHANGE HERE IF NEW GENE LIST ADDED
# -----------------------------------------------------------------------------
D25029_datasets <- list(
  top_candidate_genes_D25029_MAF0.01 = top_candidate_genes_D25029_MAF0.01,
  top_candidate_related_genes_D25029_MAF0.01 = top_candidate_related_genes_D25029_MAF0.01,
  candidate_genes_D25029_MAF0.01 = candidate_genes_D25029_MAF0.01,
  cilium_D25029_MAF0.01 = cilium_D25029_MAF0.01,
  ECM_interaction_D25029_MAF0.01 = ECM_interaction_D25029_MAF0.01,
  hsa04310_Wnt_D25029_MAF0.01 = hsa04310_Wnt_D25029_MAF0.01,
  hsa04340_hedgedog_D25029_MAF0.01 = hsa04340_hedgedog_D25029_MAF0.01,
  HP0006695_D25029_MAF0.01 = HP0006695_D25029_MAF0.01
)

# -----------------------------------------------------------------------------
# Summarize all datasets and store the results in a list
# -----------------------------------------------------------------------------
D25029_summarized_results <- lapply(D25029_datasets, summarize_by_gene_variant)
names(D25029_summarized_results) <- names(D25029_datasets)

# -----------------------------------------------------------------------------
# Example: Accessing the summarized result for "top_candidate_genes_D25029_MAF0.01"
# -----------------------------------------------------------------------------
print(D25029_summarized_results$top_candidate_genes_D25029_MAF0.01, n = Inf)
