# -------------------------------------------------------------------
# Focus on gene and investigate Gene.refgene
# -------------------------------------------------------------------

# 1. Gene.refgene within group (all genes)
venn_data_all <- list(
  DS_ECD = c(D25029_MAF0.01$Gene.refgene, D25046_MAF0.01$Gene.refgene),
  DS_nonECD = c(D25163_MAF0.01$Gene.refgene, D25165_MAF0.01$Gene.refgene, D25168_MAF0.01$Gene.refgene),
  nonDS_ECD = c(D25007_MAF0.01$Gene.refgene)
)

# 2. Gene.refgene within group (focus on exonic, UTR variants, and splicing-affecting variants)
venn_data_eu <- list(
  DS_ECD = c(D25029_variantgroup$exonic$Gene.refgene, D25029_variantgroup$UTR$Gene.refgene, D25029_GeneList$Af_splicing$Gene.refgene,
             D25046_variantgroup$exonic$Gene.refgene, D25046_variantgroup$UTR$Gene.refgene, D25046_GeneList$Af_splicing$Gene.refgene),
  
  DS_nonECD = c(D25163_variantgroup$exonic$Gene.refgene, D25163_variantgroup$UTR$Gene.refgene, D25163_GeneList$Af_splicing$Gene.refgene,
                D25165_variantgroup$exonic$Gene.refgene, D25165_variantgroup$UTR$Gene.refgene, D25165_GeneList$Af_splicing$Gene.refgene,
                D25168_variantgroup$exonic$Gene.refgene, D25168_variantgroup$UTR$Gene.refgene, D25168_GeneList$Af_splicing$Gene.refgene),
  
  nonDS_ECD = c(D25007_variantgroup$exonic$Gene.refgene, D25007_variantgroup$UTR$Gene.refgene, D25007_GeneList$Af_splicing$Gene.refgene)
)

# 3. Venn diagram for visualizing the groups
ggVennDiagram(venn_data_eu, color = "black", 
              lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_x_continuous(expand = expansion(mult = .2))

# 4. Exclusively Gene.refgene in DS-ECD, excluding those shared with DS_nonECD
exclusive_refgene_DSECDID <- setdiff(venn_data_eu$DS_ECD, venn_data_eu$DS_nonECD)

# 5. Extract the exclusive variants from the datasets (exonic, UTR)
exclusive_refgene_DSECD <- rbind(
  D25029_variantgroup$exonic %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25029_variantgroup$UTR %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25046_variantgroup$exonic %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25046_variantgroup$UTR %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID)
) %>% distinct(Gene.refgene, .keep_all = TRUE)

# 6. Extract exclusive variants from the Af_splicing group
exclusive_refgene_DSECD_AFSPL <- rbind(
  D25029_GeneList$Af_splicing %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25046_GeneList$Af_splicing %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID)
) %>% distinct(Gene.refgene, .keep_all = TRUE)

# 7. Save the results to Excel files
wb <- createWorkbook()
addWorksheet(wb, "exclusive_refgene_DSECD")
writeData(wb, "exclusive_refgene_DSECD", exclusive_refgene_DSECD, startRow = 1, startCol = 1)
saveWorkbook(wb, "exclusive_refgene_DSECD.xlsx", overwrite = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "exclusive_refgene_DSECD_AFSPL")
writeData(wb, "exclusive_refgene_DSECD_AFSPL", exclusive_refgene_DSECD_AFSPL, startRow = 1, startCol = 1)
saveWorkbook(wb, "exclusive_refgene_DSECD_AFSPL.xlsx", overwrite = TRUE)
