#### Fetch Gene List of KEGG pathway ####
## Retrieve the gene list using the KEGG REST API and process it manually
pathway <- "hsa04020" # calcium signaling pathway in humans
url <- paste("http://rest.kegg.jp/link/hsa/",pathway, sep="")
gene_data <- read.table(url, header = FALSE, sep = "\t")
# Extract the KEGG gene IDs by removing the 'hsa:' prefix
gene_ids <- sub("hsa:", "", gene_data$V2)

## Convert Gene IDs to Gene Symbols
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

gene_symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
# View the gene symbols, key is geneID and value is geneSymbol
print(gene_symbols)
# Create a data frame with one column "SYMBOL"
gene_df <- data.frame(SYMBOL = gene_symbols, stringsAsFactors = FALSE)
# Write the data frame to a text file
filename <- "hsa04020.txt"
write.table(gene_df, file = filename, row.names = FALSE, col.names = FALSE, quote = FALSE)
