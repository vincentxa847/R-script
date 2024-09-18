# Summarize the data by Gene and Variant Type
summary <- function(df){
  
  summary_table <- df %>%
    group_by(Gene_refgene, Func_refgene) %>%
    summarize(count = n(), .groups = "drop")
  
  wide_table <- summary_table %>%
    tidyr::pivot_wider(names_from = Func_refgene, values_from = count, values_fill = 0)
  
  # View the summarized table
  print(wide_table, n = Inf)
  
}

# D25007
sum_top_candidate_genes_D25007_MAF0.01 = summary(top_candidate_genes_D25007_MAF0.01)
sum_top_candidate_related_genes_D25007_MAF0.01 = summary(top_candidate_related_genes_D25007_MAF0.01)
sum_candidate_genes_D25007_MAF0.01 = summary(candidate_genes_D25007_MAF0.01)
sum_cilium_D25007_MAF0.01 = summary(cilium_D25007_MAF0.01)
sum_CL8946_folic_acid_D25007_MAF0.01 = summary(CL8946_folic_acid_D25007_MAF0.01)
sum_ECM_interaction_D25007_MAF0.01 = summary(ECM_interaction_D25007_MAF0.01)
sum_hsa04310_Wnt_D25007_MAF0.01 = summary(hsa04310_Wnt_D25007_MAF0.01)
sum_hsa04340_hedgedog_D25007_MAF0.01 = summary(hsa04340_hedgedog_D25007_MAF0.01)
sum_HP0006695_D25007_MAF0.01 = summary(HP0006695_D25007_MAF0.01)

# N1675
sum_top_candidate_genes_N1675_MAF0.01 = summary(top_candidate_genes_N1675_MAF0.01)
sum_top_candidate_related_genes_N1675_MAF0.01 = summary(top_candidate_related_genes_N1675_MAF0.01)
sum_candidate_genes_N1675_MAF0.01 = summary(candidate_genes_N1675_MAF0.01)
sum_cilium_N1675_MAF0.01 = summary(cilium_N1675_MAF0.01)
sum_CL8946_folic_acid_N1675_MAF0.01 = summary(CL8946_folic_acid_N1675_MAF0.01)
sum_ECM_interaction_N1675_MAF0.01 = summary(ECM_interaction_N1675_MAF0.01)
sum_hsa04310_Wnt_N1675_MAF0.01 = summary(hsa04310_Wnt_N1675_MAF0.01)
asum_hsa04340_hedgedog_N1675_MAF0.01 = summary(hsa04340_hedgedog_N1675_MAF0.01)
sum_HP0006695_N1675_MAF0.01 = summary(HP0006695_N1675_MAF0.01)

# AI3008
sum_top_candidate_genes_AI3008_MAF0.01 = summary(top_candidate_genes_AI3008_MAF0.01)
sum_top_candidate_related_genes_AI3008_MAF0.01 = summary(top_candidate_related_genes_AI3008_MAF0.01)
sum_candidate_genes_AI3008_MAF0.01 = summary(candidate_genes_AI3008_MAF0.01)
sum_cilium_AI3008_MAF0.01 = summary(cilium_AI3008_MAF0.01)
sum_CL8946_folic_acid_AI3008_MAF0.01 = summary(CL8946_folic_acid_AI3008_MAF0.01)
sum_ECM_interaction_AI3008_MAF0.01 = summary(ECM_interaction_AI3008_MAF0.01)
sum_hsa04310_Wnt_AI3008_MAF0.01 = summary(hsa04310_Wnt_AI3008_MAF0.01)
sum_hsa04340_hedgedog_AI3008_MAF0.01 = summary(hsa04340_hedgedog_AI3008_MAF0.01)
sum_HP0006695_AI3008_MAF0.01 = summary(HP0006695_AI3008_MAF0.01)
