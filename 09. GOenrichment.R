# -----------------------------------------------------------------------------
# Using clusterProfiler to perform Over-representation analysis (ORA)
# Reference : Xu, S., Hu, E., Cai, Y. et al. Using clusterProfiler to characterize multiomics data. Nat Protoc (2024).
# -----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#### Function to perform ORA of KEGG or GO and visualize the result ####
ORA.Enrichment <- function(geneset, analysisType = "KEGG", ont = "BP", pvalueCutoff = 0.05, showCategory = 10, showCategoryNetwork = 5) {
  # Sample gene list
  gene_list <- geneset
  
  # Convert gene symbols to Entrez IDs
  gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_ids <- gene_entrez$ENTREZID
  
  # Perform KEGG or GO enrichment analysis
  if (analysisType == "KEGG") {
    # KEGG Enrichment Analysis
    enrichment_result <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = pvalueCutoff)
  } else if (analysisType == "GO") {
    # GO Enrichment Analysis (for BP, MF, or CC)
    enrichment_result <- enrichGO(gene = entrez_ids,
                                  OrgDb = org.Hs.eg.db,
                                  ont = ont, # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                                  pvalueCutoff = pvalueCutoff)
  } else {
    stop("Invalid analysis type. Choose either 'KEGG' or 'GO'.")
  }
  
  # Create barplot and save as object
  bar_plot <- barplot(enrichment_result, showCategory = showCategory)
  
  # Create dotplot and save as object
  dot_plot <- dotplot(enrichment_result, showCategory = showCategory)
  
  # Network plot: Convert Entrez IDs back to gene symbols
  gene_symbol_df <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  
  # Replace Entrez IDs with Gene Symbols in enrichment results
  enrichment_result@result$geneID <- sapply(enrichment_result@result$geneID, function(x) {
    ids <- unlist(strsplit(x, "/"))
    symbols <- gene_symbol_df$SYMBOL[match(ids, gene_symbol_df$ENTREZID)]
    paste(symbols, collapse = "/")
  })
  
  # Create cnetplot and save as object
  network_plot <- cnetplot(enrichment_result, showCategory = showCategoryNetwork)
  
  # Return a list with enrichment results and all plots
  return(list(
    enrichment = enrichment_result,
    barplot = bar_plot,
    dotplot = dot_plot,
    cnetplot = network_plot
  ))
}

# Example usage for KEGG enrichment
tmp_kegg <- ORA.Enrichment(D25046_variantgroup$intronic$Gene_refgene, analysisType = "KEGG", pvalueCutoff = 0.05, showCategory = 10, showCategoryNetwork = 5)

# Example usage for GO enrichment (Biological Process)
tmp_go <- ORA.Enrichment(D25046_variantgroup$intronic$Gene_refgene, analysisType = "GO", ont = "BP", pvalueCutoff = 0.05, showCategory = 10, showCategoryNetwork = 2)

# Access the barplot using $ operator
tmp_kegg$barplot
tmp_go$barplot

# Access the dotplot using $ operator
tmp_kegg$dotplot
tmp_go$dotplot

# Access the cnetplot using $ operator
tmp_kegg$cnetplot
tmp_go$cnetplot

# View the enrichment results
head(tmp_kegg$enrichment)
head(tmp_go$enrichment)
