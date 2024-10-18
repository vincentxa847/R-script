# Load necessary database
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

# Function to fetch gene list from KEGG pathway
fetch_kegg_gene_list <- function(pathway_id, output_file) {
  # Construct the URL to fetch gene data from KEGG API
  # TODO: If http cannot fetch then fall back to https here
  url <- paste0("https://rest.kegg.jp/link/hsa/", pathway_id)
  
  # Retrieve gene data from the URL
  gene_data_raw <- readLines(url)
  gene_data <- read.table(text = paste(gene_data_raw, collapse = "\n"), header = FALSE, sep = "\t")
  
  # Extract the KEGG gene IDs by removing the 'hsa:' prefix
  gene_ids <- sub("hsa:", "", gene_data$V2)
  
  # Convert Gene IDs to Gene Symbols
  gene_symbols <- mapIds(
    org.Hs.eg.db, 
    keys = gene_ids, 
    column = "SYMBOL", 
    keytype = "ENTREZID", 
    multiVals = "first"
  )
  
  # Create a data frame with one column "SYMBOL"
  gene_df <- data.frame(SYMBOL = gene_symbols, stringsAsFactors = FALSE)
  
  # Write the data frame to a text file
  write.table(
    gene_df, 
    file = output_file, 
    row.names = FALSE, 
    col.names = FALSE, 
    quote = FALSE
  )
  
  # Return the gene data frame for inspection
  return(gene_df)
}

# Usage of the function
pathway_id <- "hsa04020" # Calcium signaling pathway in humans
output_file <- "hsa04020.txt"

# Fetch the gene list and save it to a file
gene_list_df <- fetch_kegg_gene_list(pathway_id, output_file)

# Print the gene symbols
print(gene_list_df)
