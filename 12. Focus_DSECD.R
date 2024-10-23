#### 1. Gene.refgene within group (all genes)
venn_data_all <- list(
  DS_ECD = c(D25029_MAF0.01$Gene.refgene, D25046_MAF0.01$Gene.refgene),
  DS_nonECD = c(D25163_MAF0.01$Gene.refgene, D25165_MAF0.01$Gene.refgene, D25168_MAF0.01$Gene.refgene),
  nonDS_ECD = c(D25007_MAF0.01$Gene.refgene)
)

#### 2. Gene.refgene within group (focus on exonic, UTR variants, and splicing-affecting variants)
venn_data_eu <- list(
  DS_ECD = c(D25029_variantgroup$exonic$Gene.refgene, D25029_variantgroup$UTR$Gene.refgene, D25029_GeneList$Af_splicing$Gene.refgene,
             D25046_variantgroup$exonic$Gene.refgene, D25046_variantgroup$UTR$Gene.refgene, D25046_GeneList$Af_splicing$Gene.refgene),
  
  DS_nonECD = c(D25163_variantgroup$exonic$Gene.refgene, D25163_variantgroup$UTR$Gene.refgene, D25163_GeneList$Af_splicing$Gene.refgene,
                D25165_variantgroup$exonic$Gene.refgene, D25165_variantgroup$UTR$Gene.refgene, D25165_GeneList$Af_splicing$Gene.refgene,
                D25168_variantgroup$exonic$Gene.refgene, D25168_variantgroup$UTR$Gene.refgene, D25168_GeneList$Af_splicing$Gene.refgene),
  
  nonDS_ECD = c(D25007_variantgroup$exonic$Gene.refgene, D25007_variantgroup$UTR$Gene.refgene, D25007_GeneList$Af_splicing$Gene.refgene)
)

#### 3. Venn diagram for visualizing the groups
ggVennDiagram(venn_data_eu, color = "black", 
              lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
  scale_x_continuous(expand = expansion(mult = .2))

#### 4. Exclusively Gene.refgene in DS-ECD, excluding those shared with DS_nonECD
exclusive_refgene_DSECDID <- setdiff(venn_data_eu$DS_ECD, venn_data_eu$DS_nonECD)

# Intersect of DS_ECD and nonDS_ECD
shared_genes_DSECD_nonDSECD <- intersect(venn_data_eu$DS_ECD, venn_data_eu$nonDS_ECD)
shared_genes_all_three <- Reduce(intersect, venn_data_eu)
exclusive_shared_DSECD_nonDSECD <- setdiff(shared_genes_DSECD_nonDSECD, shared_genes_all_three)


#### 5. Extract the variants (exonic, UTR, splicing) in exclusive genes from the datasets 
exclusive_refgene_DSECD <- rbind(
  D25029_variantgroup$exonic[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25029_variantgroup$UTR[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25029_GeneList$Af_splicing[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25046_variantgroup$exonic[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25046_variantgroup$UTR[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID),
  D25046_GeneList$Af_splicing[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_refgene_DSECDID)
) %>% distinct(Gene.refgene, .keep_all = TRUE)

exclusive_shared_DSECD_nonDSECD <- rbind(
  D25029_variantgroup$exonic[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_shared_DSECD_nonDSECD),
  D25029_variantgroup$UTR[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_shared_DSECD_nonDSECD),
  D25029_GeneList$Af_splicing[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_shared_DSECD_nonDSECD),
  D25046_variantgroup$exonic[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_shared_DSECD_nonDSECD),
  D25046_variantgroup$UTR[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_shared_DSECD_nonDSECD),
  D25046_GeneList$Af_splicing[,c(1:27)] %>% filter(Gene.refgene %in% exclusive_shared_DSECD_nonDSECD)
) %>% distinct(Gene.refgene, .keep_all = TRUE)


#### 6. Save the results to Excel files
wb <- createWorkbook()
addWorksheet(wb, "exclusivegene_DSECD")
writeData(wb, "exclusivegene_DSECD", exclusive_refgene_DSECD, startRow = 1, startCol = 1)
saveWorkbook(wb, "../EXCLUDE_DS_GENE/exclusivegene_DSECD.xlsx", overwrite = TRUE)

wb <- createWorkbook()
addWorksheet(wb, "exclusive_shared_DSECD_nonDSECD")
writeData(wb, "exclusive_shared_DSECD_nonDSECD", exclusive_shared_DSECD_nonDSECD, startRow = 1, startCol = 1)
saveWorkbook(wb, "../EXCLUDE_DS_GENE/exclusive_shared_DSECD_nonDSECD.xlsx", overwrite = TRUE)

