library(org.Hs.eg.db)     

# -----------------------------------------------------------------------------
# Plot genes associated with variants in the KEGG pathway. 
# This visualization includes all genes with variant(s), regardless of the type or number of variants they have.
# Genes with available CADD scores will be highlighted in the plot.
# -----------------------------------------------------------------------------
kegg_genesWithVariant_CADD <- function(sample, table_name, pathway_id, output_dir) {
  # Extract relevant columns: Gene_refgene and CADD_phred
  gene_cadd_df <- sample[[table_name]][, c("Gene_refgene", "CADD_phred")]
  # Handle CADD value
  gene_cadd_df$CADD_phred[gene_cadd_df$CADD_phred == "."] <- 0
  gene_cadd_df$CADD_phred <- as.numeric(gene_cadd_df$CADD_phred)
  
  # Remove duplicates based on Gene_refgene and CADD_phred
  unique_gene_cadd_df <- unique(gene_cadd_df)
  
  # Convert gene symbols to Entrez IDs
  entrez_genes <- clusterProfiler::bitr(unique_gene_cadd_df$Gene_refgene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Merge Entrez IDs with CADD_phred values
  entrez_genes <- merge(entrez_genes, unique_gene_cadd_df, by.x = "SYMBOL", by.y = "Gene_refgene")
  
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
kegg_genesWithVariant_CADD(D25165_GeneList, "hsa04310_Wnt", "hsa04310", "../DS-nonECD/D25165/")
