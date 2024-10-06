library(dplyr)

#### Prepare gene list ####
# Function to read gene lists from text files and handle errors
read_gene_list <- function(file_path) {
  if (file.exists(file_path)) {
    if (grepl("HPOList", file_path)) {
      # For files in "HPOList", read with header and select "name" column (this is directly download from HPO website)
      gene_list <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)$name
      gene_list <- trimws(gene_list)
    } else {
      # For other files, read without a header and select the first column
      gene_list <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
    }
    return(gene_list)  
  } else {
    stop(paste("File not found:", file_path))  
  }
}

# Define paths for the gene lists
# TODO: REMAIN TO UPDATE
gene_paths <- list(
  top_candidate_genes = "../GENE_LIST/TOP_candidate_1006.txt",
  candidate_genes = "../GENE_LIST/candidate_1006.txt",
  top_candidate_related_genes = "../GENE_LIST/TOP_candidate_STRING_1006.txt",
  # GO
  cilium = "../GENE_LIST/GO/GO0005929_cilium.txt", # CHANGE TO SCGS_v2
  protein_folding = "../GENE_LIST/GO/GO0140662_ATPdependent protein folding chaperone.txt",
  # KEGG
  ECM_interaction = "../GENE_LIST/KEGG/hsa04512_ECM_receptor_interaction.txt",
  hsa04310_Wnt = "../GENE_LIST/KEGG/hsa04310_Wnt.txt",
  hsa04330_Notch = "../GENE_LIST/KEGG/hsa04330_Notch.txt",
  hsa04340_hedgedog = "../GENE_LIST/KEGG/hsa04340_hedgedog.txt",
  hsa04020_calcium_signaling = "../GENE_LIST/KEGG/hsa04020_calcium_signaling.txt",
  hsa04350_TGFbata = "../GENE_LIST/KEGG/hsa04350_TGFbata.txt",
  # HPO
  HP0006695_ECD = "../GENE_LIST/HPOList/genes_for_HP_0006695"
)

# Read in gene list
gene_lists <- lapply(gene_paths, read_gene_list)

#### Output a list of gene for IPA ####
# -----------------------------------------------------------------------------
# TODO: Investigate how to use IPA (Ingenuity Pathway Analysis) for variant analysis
#       and determine what additional information should be included.
# -----------------------------------------------------------------------------

# Function to extract gene symbols from variants and save them to a text file
# If more than one variant exists in a gene, multiple gene symbols may be written.
# These duplicates will be handled (removed) during IPA processing.
writeout_list <- function(table, gene_list, output_name) {
  output <- table %>% filter(Gene_refgene %in% gene_list)
  output <- output$Gene_refgene 
  
  file_path <- paste0("../IPA/", output_name, ".txt")
  
  write.table(output, file = file_path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Example usage : 
out <- writeout_list(D25029,gene_lists$HP0006695_ECD, "D25029_MAF0.01_HP0006695")  

#### Analyze in current section ####
# Function to see if variants in genes of genelist
gene_list_filter <- function(input, gene_list) {
  if (!"Gene_refgene" %in% colnames(input)) {
    stop("The input data does not contain a column named 'Gene_refgene'.")
  }
  
  return(input %>% filter(Gene_refgene %in% gene_list))
}

# Example usage : D25029 (DS-ECD)
top_candidate_genes_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$top_candidate_genes)
top_candidate_related_genes_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$top_candidate_related_genes)
candidate_genes_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$candidate_genes)
ECM_interaction_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$ECM_interaction)
hsa04310_Wnt_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$hsa04310_Wnt)
cilium_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$cilium)
hsa04330_Notch_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$hsa04330_Notch)
hsa04310_Wnt_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$hsa04310_Wnt)
hsa04340_hedgedog_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$hsa04340_hedgedog)
CL8946_folic_acid_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$CL8946_folic_acid)
HSA189451_heme_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$HSA189451_heme)
HP0006695_D25029_MAF0.01 = gene_list_filter(D25029,gene_lists$HP0006695)

