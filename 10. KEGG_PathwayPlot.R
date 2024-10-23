library(org.Hs.eg.db)
library(pathview)
# -----------------------------------------------------------------------------
# Plot genes associated with variants in the KEGG pathway. 
# This visualization includes all genes with variant(s), regardless of the type or number of variants they have.
# Genes with available CADD scores will be highlighted in the plot.
# -----------------------------------------------------------------------------
kegg_genesWithVariant_CADD <- function(sample, table_name, pathway_id, output_dir, isTable=FALSE) {
  
  # WHERE TABLE RATHER THAN LIST INPUT
  if (isTable) {
    # Extract relevant columns: Gene.refgene and CADD_phred
    gene_cadd_df <- sample[, c("Gene.refgene", "CADD_phred")]
  }
  else{
    # Extract relevant columns: Gene.refgene and CADD_phred
    gene_cadd_df <- sample[[table_name]][, c("Gene.refgene", "CADD_phred")]
  }
  # Handle CADD value
  gene_cadd_df$CADD_phred[gene_cadd_df$CADD_phred == "."] <- 0
  gene_cadd_df$CADD_phred <- as.numeric(gene_cadd_df$CADD_phred)
  
  # Check if all CADD_phred values are zero (no CADD score), which will cause pathview fail
  if (all(gene_cadd_df$CADD_phred == 0)) {
    message("All CADD_phred values are zero. Assigning 0.001 to the first gene for visualization.")
    gene_cadd_df$CADD_phred[1] <- 0.001  # Assign 0.001 to the first row
  }
  
  # Remove duplicates based on Gene.refgene and CADD_phred
  unique_gene_cadd_df <- unique(gene_cadd_df)
  
  # Convert gene symbols to Entrez IDs
  entrez_genes <- clusterProfiler::bitr(unique_gene_cadd_df$Gene.refgene, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db)
  
  # Merge Entrez IDs with CADD_phred values
  entrez_genes <- merge(entrez_genes, unique_gene_cadd_df, by.x = "SYMBOL", by.y = "Gene.refgene")
  
  # Filter rows where CADD_phred is not zero
  non_zero_cadd_genes <- entrez_genes[entrez_genes$CADD_phred != 0, ]
  
  # Print the SYMBOL, ENTREZID, and CADD_phred values
  message("Genes with non-zero CADD_phred values:")
  print(non_zero_cadd_genes[, c("SYMBOL", "ENTREZID", "CADD_phred")])
  
  # Create a named vector with Entrez IDs as names and CADD_phred as values
  gene_data <- setNames(entrez_genes$CADD_phred, entrez_genes$ENTREZID)
  
  # Store the current working directory
  original_wd <- getwd()
  
  # Change the working directory to the specified output directory
  output_dir <- paste0(output_dir, pathway_id, ".pathview")
  dir.create(output_dir, recursive = FALSE, showWarnings = FALSE)  # Create directory if it doesn't exist
  setwd(output_dir)  # Set the new working directory
  
  # Visualize the pathway, using CADD_phred scores for highlighting the genes
  pathview::pathview(gene.data = gene_data, pathway.id = pathway_id, species = "hsa", gene.idtype = "ENTREZID", 
                     limit = list(gene = c(min(gene_data, na.rm = TRUE), max(gene_data, na.rm = TRUE))))
  
  # Revert back to the original working directory
  setwd(original_wd)
}

# -----------------------------------------------------------------------------
#### Usage ####
# -----------------------------------------------------------------------------
## nonDS-ECD
# D25007
kegg_genesWithVariant_CADD(D25007_GeneList, "hsa04310_Wnt", "hsa04310", "../nonDS-ECD/D25007/")
kegg_genesWithVariant_CADD(D25007_GeneList, "hsa04330_Notch", "hsa04330", "../nonDS-ECD/D25007/")
kegg_genesWithVariant_CADD(D25007_GeneList, "hsa04340_hedgedog", "hsa04340", "../nonDS-ECD/D25007/")
kegg_genesWithVariant_CADD(D25007_GeneList, "hsa04020_calcium_signaling", "hsa04020", "../nonDS-ECD/D25007/")
kegg_genesWithVariant_CADD(D25007_GeneList, "hsa04350_TGFbeta", "hsa04350", "../nonDS-ECD/D25007/")
#kegg_genesWithVariant_CADD(D25007_GeneList, "N01453_BMP", "N01453", "../nonDS-ECD/D25007/")
## DS-nonECD
# D25165
kegg_genesWithVariant_CADD(D25165_GeneList, "hsa04310_Wnt", "hsa04310", "../DS-nonECD/D25165/")
kegg_genesWithVariant_CADD(D25165_GeneList, "hsa04330_Notch", "hsa04330", "../DS-nonECD/D25165/")
kegg_genesWithVariant_CADD(D25165_GeneList, "hsa04340_hedgedog", "hsa04340", "../DS-nonECD/D25165/")
kegg_genesWithVariant_CADD(D25165_GeneList, "hsa04020_calcium_signaling", "hsa04020", "../DS-nonECD/D25165/")
kegg_genesWithVariant_CADD(D25165_GeneList, "hsa04350_TGFbeta", "hsa04350", "../DS-nonECD/D25165/")
#kegg_genesWithVariant_CADD(D25165_GeneList, "N01453_BMP", "N01453", "../DS-nonECD/D25165/")
# D25168
kegg_genesWithVariant_CADD(D25168_GeneList, "hsa04310_Wnt", "hsa04310", "../DS-nonECD/D25168/")
kegg_genesWithVariant_CADD(D25168_GeneList, "hsa04330_Notch", "hsa04330", "../DS-nonECD/D25168/")
kegg_genesWithVariant_CADD(D25168_GeneList, "hsa04340_hedgedog", "hsa04340", "../DS-nonECD/D25168/")
kegg_genesWithVariant_CADD(D25168_GeneList, "hsa04020_calcium_signaling", "hsa04020", "../DS-nonECD/D25168/")
kegg_genesWithVariant_CADD(D25168_GeneList, "hsa04350_TGFbeta", "hsa04350", "../DS-nonECD/D25168/")
## DS-ECD
# D25029 
kegg_genesWithVariant_CADD(D25029_GeneList, "hsa04310_Wnt", "hsa04310", "../DS-ECD/D25029/")
kegg_genesWithVariant_CADD(D25029_GeneList, "hsa04330_Notch", "hsa04330", "../DS-ECD/D25029/")
kegg_genesWithVariant_CADD(D25029_GeneList, "hsa04340_hedgedog", "hsa04340", "../DS-ECD/D25029/")
kegg_genesWithVariant_CADD(D25029_GeneList, "hsa04020_calcium_signaling", "hsa04020", "../DS-ECD/D25029/")
kegg_genesWithVariant_CADD(D25029_GeneList, "hsa04350_TGFbeta", "hsa04350", "../DS-ECD/D25029/")
# D25046 
kegg_genesWithVariant_CADD(D25046_GeneList, "hsa04310_Wnt", "hsa04310", "../DS-ECD/D25046/")
kegg_genesWithVariant_CADD(D25046_GeneList, "hsa04330_Notch", "hsa04330", "../DS-ECD/D25046/")
kegg_genesWithVariant_CADD(D25046_GeneList, "hsa04340_hedgedog", "hsa04340", "../DS-ECD/D25046/") 
kegg_genesWithVariant_CADD(D25046_GeneList, "hsa04020_calcium_signaling", "hsa04020", "../DS-ECD/D25046/")
kegg_genesWithVariant_CADD(D25046_GeneList, "hsa04350_TGFbeta", "hsa04350", "../DS-ECD/D25046/")

# Exclusively Gene.refgene in DS-ECD and DS-nonECD
kegg_genesWithVariant_CADD(exclusive_refgene_DSECD, "hsa04310_Wnt", "hsa04310", "../EXCLUDE_DS_GENE/exclusive_DSECD_NOTDSnonECD/",TRUE)
kegg_genesWithVariant_CADD(exclusive_refgene_DSnonECD, "hsa04310_Wnt", "hsa04310", "../EXCLUDE_DS_GENE/exclusive_DSnonECD_NOTDSECD/",TRUE)

kegg_genesWithVariant_CADD(exclusive_refgene_DSECD, "hsa04330_Notch", "hsa04330", "../EXCLUDE_DS_GENE/exclusive_DSECD_NOTDSnonECD/",TRUE)
kegg_genesWithVariant_CADD(exclusive_refgene_DSnonECD, "hsa04330_Notch", "hsa04330", "../EXCLUDE_DS_GENE/exclusive_DSnonECD_NOTDSECD/",TRUE)

kegg_genesWithVariant_CADD(exclusive_refgene_DSECD, "hsa04340_hedgedog", "hsa04340", "../EXCLUDE_DS_GENE/exclusive_DSECD_NOTDSnonECD/",TRUE)
kegg_genesWithVariant_CADD(exclusive_refgene_DSnonECD, "hsa04340_hedgedog", "hsa04340", "../EXCLUDE_DS_GENE/exclusive_DSnonECD_NOTDSECD/",TRUE)

kegg_genesWithVariant_CADD(exclusive_refgene_DSECD, "hsa04020_calcium_signaling", "hsa04020", "../EXCLUDE_DS_GENE/exclusive_DSECD_NOTDSnonECD/",TRUE)
kegg_genesWithVariant_CADD(exclusive_refgene_DSnonECD, "hsa04020_calcium_signaling", "hsa04020", "../EXCLUDE_DS_GENE/exclusive_DSnonECD_NOTDSECD/",TRUE)

kegg_genesWithVariant_CADD(exclusive_refgene_DSECD, "hsa04350_TGFbeta", "hsa04350", "../EXCLUDE_DS_GENE/exclusive_DSECD_NOTDSnonECD/",TRUE)
kegg_genesWithVariant_CADD(exclusive_refgene_DSnonECD, "hsa04350_TGFbeta", "hsa04350", "../EXCLUDE_DS_GENE/exclusive_DSnonECD_NOTDSECD/",TRUE)

kegg_genesWithVariant_CADD(exclusive_refgene_DSECD, "hsa04390_Hippo", "hsa04390", "../EXCLUDE_DS_GENE/exclusive_DSECD_NOTDSnonECD/",TRUE)
kegg_genesWithVariant_CADD(exclusive_refgene_DSnonECD, "hsa04390_Hippo", "hsa04390", "../EXCLUDE_DS_GENE/exclusive_DSnonECD_NOTDSECD/",TRUE)
